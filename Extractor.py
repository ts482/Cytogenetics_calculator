#importing modules

import re
import pandas as pd
#import streamlit as st
#import base64
punct = r'()[]/'
punct_dict = {s:0 for s in punct}


def load_file(file):
    '''
    Loads the excel file and drops the ID column
    
    Params:
    -------
    
    file: excel filepath
        the file containing the cytogenetic reports, with columns
        that contain information about details to fill
        
    streamlit: bool, default False
        whether the streamlit interface should be used
    '''

    karyotypes = pd.read_csv(file)    
    karyotypes = karyotypes.drop(columns='ID')
    
    return karyotypes

def properties_dict(karyotypes, properties = None):
    '''
    1. Transforms column names in input file into relevant 
    string formats to be used in regex report extraction
    2. Returns a dictionary having these strings as keys,
    with the original column names as values
    
    Params:
    -------
    karyotypes: pd.DataFrame
        dataframe obtained in load_file
        
    properties: list, default None
        option to obtain string names for abnormalities directly.
    '''
    if not properties:
        properties = karyotypes.columns[4:-1].to_list()
    d = {p:p for p in properties}
    for k,v in d.items():
        #removing quotation marks at start and end
        if v[0] == '"':
            if v[-1] == v[0]:
                v = v[1:-1]
        
        #Formatting monosomies
        if v.startswith('Monosomy'):
            v = '-' + v[8:]
        
        #
        matches = re.finditer('(\d+)(p|q)',v)
        for m in matches:
            #different rule for translocations - removing first bracket
            if re.search('t\(', m.string):
                f = v[:m.start()] + m.group()[:-1] + \
                    ')(' + m.group()[-1]+ v[m.end():]
            else:
                f = v[:m.start()] + '(' + m.group()[:-1] + \
                    ')(' + m.group()[-1]+ v[m.end():]
            v=f
        
        
            
        #creating escape characters for strings
        v = re.escape(v)
        
        #special case for dots
        if v == "t\(8;16\)\(p11;p13\)":
            v = "t\(8;16\)\(p11\.?\d?;p13\.?\d?\)"
        
        #creating special case for for t(v;11)
        if v == 't\\(v;11\\)':
            v = '(t\\(\d+;11\\))|(t\\(11;\d+\\))'
        #another special case for t(3q26.2;v)
        if v == 't\\(3\\)\\(q26\\.2;v\\)':
            v = '(t\(\d+;3\)\([pq]\d+(\.\d+)?;q26\.2\))|(t\(3;\d+\)\(q26\.2;[pq]\d+(\.\d+)?\))'

        
        d[k] = v
    return {v:k for k,v in d.items()}


def remove_artefact(row):
    '''
    Removes whitespaces and import artefacts from reports 
    Function to be used in conjunction with pd.DataFrame.apply, 
    and the column containing the cytogenetic report
    '''
    if pd.notna(row['Cytogenetics']):
        if re.search('_x000D_$', row['Cytogenetics']):
            return row['Cytogenetics'].replace(' ', '')[:-7]
        return row['Cytogenetics'].replace(' ', '')
    return 'Error'


def gram_error(string, verbose=False):
    '''
    Takes in a cytogenetic report string and checks for multiple error types.
    Returns a list describing the errors picked up.
    
    Error types included:
    1. Automated process failure (e.g. presence of 'Error' or 'fail')
    2. Warnings in presence of constitutional changes
    3. Deducible punctuation rules, usually to do with ratios of sections
    4. Counting expected number of chromosomes and comparing to stated number
    for each report section
    '''
    error = []
    warning = []
    chr_count = {}
    
    #basic report errors
    if string == 'Error':
        error.append('String report missing')
        
    if re.search('fail', string.lower()):
        error.append('String report indicates failure')
    
    #highlight if error contains uncertainty
    if re.search('\\?', string):
        warning.append('Question mark in string means report is uncertain')
        
    #all reports should include commas somewhere 
    if not re.search(',', string):
        error.append('Missing comma')

    #checking for YX

    if re.search('YX', string):
        warning.appened('"YX" appears to be wrong way round. Did you mean "XY"?')
        
    #checking for constitutional changes
    if re.search('[^a-z]?c[^p]', string) or re.search('[^a-z]c$', string):
        warning.append('constitutional changes present')
        
    #incorrectly counts abnormalities if split character are together
    if re.search('[/,[]{2}', string):
        error.append('The following symbols must not be next to another or themselves: / , [')
    
    #grammatical rules regarding symbol ratios
    missing = list()
    for k in punct_dict:
        punct_dict[k] = len(re.findall(re.escape(k), string))
    if punct_dict['/'] != punct_dict['[']-1:
        error.append('incorrrect ratio of \'\\\' to \'[]\' ')
    if punct_dict['['] != punct_dict[']']:
        missing.append(min(['[', ']'], key=punct_dict.get))
    if punct_dict['('] != punct_dict[')']:
        missing.append(min(['(',')'], key=punct_dict.get))
    if len(missing)!=0:
        error.append(f'missing grammar: {", ".join(missing)}')
        
        
    #checking number in subsections matches stated number
    substring = re.split('/', string)
    for i, s in enumerate(substring):
        try:
            chrom = s[:s.index(',')] #chromosome string
        except ValueError:
            error.append('Part of report missing comma')
            continue
        if not re.search('idem', s):
            expected = 46
        expected -= len(re.findall('\-', s))
        expected += len(re.findall('\+', s))
        if re.search('mar', s):
            mar_plural = re.search('\+(\d)(~\d)?mar',s)
            if mar_plural:
                expected += int(mar_plural.groups()[0]) - 1
                tilda_present =  re.search('\+(\d)~(\d)mar', s)
                #error if variable number of markers included (e.g. +2~4mar)
                if tilda_present:
                    d = tilda_present.group(1)
                    warning.append(f'Variable number of markers in report. Minimum number ({d}) used by default')
        if verbose:
            chr_count['Expected number'] = expected
        
        #in the case of a range of potential chromosome counts
        if re.search('[~]', chrom):
            try:
                low_num, high_num = int(chrom[:2]), int(chrom[-2:])
            except ValueError:
                error.append('Part of report not clearly defined by two chromosome numbers followed by comma (e.g. "43~45,")')
                continue
            if verbose:
                chr_count['Reported number'] = str(low_num) + '~' + str(high_num)
            if low_num <= expected <= high_num:
                pass
            elif expected > high_num:
                warning.append(f'chromosome number lower than expected in subsection number {i+1}')
            else:
                warning.append(f'chromosome number higher than expected in subsection number {i+1}')
                
        #if only one chromosome count is present
        else:
            try:
                num = int(chrom[:2])
            except ValueError:
                error.append('Start of report missing clear chromosome number followed by comma (e.g. "46,")')
                continue
            if verbose:
                chr_count['Reported number'] = num
            if expected == num:
                pass
            elif expected < num:
                warning.append(f'chromosome number higher than expected in subsection number {i+1}')
            else:
                warning.append(f'chromosome number lower than expected in subsection number {i+1}')
    return error, warning, chr_count

def make_multi_translocation_dict(string):
    '''
    Creates a translocation list if only full chromosomes are involved
    in translocation, or a dictionary if chromosome arms are involved
    '''
    
    #finding out how many chromosomes in translocation
    exp = '(\d\d?);(\d\d?)'
    while re.search(exp, string):
        chr_groups = list(re.search(exp, string).groups())
        exp = exp + ';(\d\d?)'
        
    #adding arms if present
    if re.search('[pq]', string):
        exp2 = '([pq]\d\d?\.?\d?\d?)'
        for g in range(len(chr_groups)-1):
            exp2 = exp2 + ';([pq]\d\d?\.?\d?\d?)'
        arm_groups = list(re.search(exp2, string).groups())
        return {int(c):a for c,a in zip(chr_groups, arm_groups)}
    return chr_groups

def check_trans_dict(trans_dict, col_true, verbose=False):
    '''
    Takes in multi translocation dictionary and detects abnormalities associated
    with adjacent chromosome translocation
    '''
    if verbose:
        verbose_list = []
    chr_keys = list(trans_dict)
    for i in range(len(chr_keys)-1):
        sorted_list = sorted([chr_keys[i],chr_keys[i+1]])
        first_chr, second_chr = sorted_list[0], sorted_list[1]
        for p in prop_dict:
            if re.search(p, f't({first_chr};{second_chr})'):
                col_true.add(p)
                if verbose:
                    verbose_list.append(p)
    if verbose:
        return col_true, verbose_list
    return col_true

def two_part_translocation_check(two_part_translocation,col_true, prop_dict,
                                 verbose=False):
    if verbose:
        verbose_list = []
    match = re.search('t\(([XY\d]{1,2});([XY\d]{1,2})\)\(([pq])\d*;([pq])\d*\)',
              two_part_translocation)
    first_string = 't(' + match.group(1) + ')(' + match.group(3) + ')'
    second_string = 't(' + match.group(2) + ')(' + match.group(4) + ')'
    for p in prop_dict:
        if re.search(p, first_string):
            col_true.add(p)
            if verbose:
                verbose_list.append(p)
        if re.search(p, second_string):
            col_true.add(p)
            if verbose:
                verbose_list.append(p)
    if verbose:
        return col_true, verbose_list
    return col_true
    

def parse_karyotype(row, prop_dict, verbose=False):
    '''
    Working row by row on a dataframe, detects abnormalities
    and fills in boolean and integer counts.
    '''
    error, warning, chr_count = gram_error(row['Cytogenetics'])
    if verbose:
        row['chr_count'] = chr_count
        
    #including warnings
    row['Warnings'] = warning
    
    #checking for error
    if error:
        row['Error description'] = error
        row['Error'] = True
        return row
    
    verbose_dict = {}

    #preparing string before splitting to abnormalities
    full_string = row['Cytogenetics']
    full_string = re.sub('\?', '', full_string)
    
    #initialising counts and sets for abnormalities
    abnorms = set(re.split('/|,|\[', full_string))
    removed = set()
    col_true = set()
    mono = 0
    poly = 0
    struc = 0
    der = 0
    mar = 0
    seventeen_p = False
    for a in abnorms:
        if verbose:
            verbose_dict[a] = []
        #cremoving normal report segments
        if re.fullmatch('(\d\d([~-]\d\d)?|[XY][XY]?|(cp)?\d\d?\]|idem)(\??c)?', a):
            removed.add(a)
            if verbose:
                if a.endswith(']'):
                    del verbose_dict[a]
                else:
                    verbose_dict[a] = 'Normal'
#         if re.search('mar', a): #no longer the case to remove mar
#             removed.add(a)

        #checking abnormalities
        if a not in removed:
            
            #checking multi-translocations
            if re.search('t\(\d\d?;\d\d?;\d\d?',a):
                trans_dict = make_multi_translocation_dict(a)
                if verbose:
                    col_true, verbose_list = \
                        check_trans_dict(trans_dict, col_true, verbose=True)
                    for v in verbose_list:
                        verbose_dict[a].append(v)
                else:
                    col_true = check_trans_dict(trans_dict, col_true)

            #checking two-part translocations
            if re.search('t\(([XY\d]{1,2});([XY\d]{1,2})\)\(([pq])\d*;([pq])\d*\)',
                         a):
                if verbose:
                    col_true, verbose_list = two_part_translocation_check(a,
                        col_true, prop_dict, verbose=True)
                    for v in verbose_list:
                        verbose_dict[a].append(v)
                else:
                    col_true = two_part_translocation_check(a,
                                col_true, prop_dict, verbose-False)
                
            
            #counting number of markers present
            if re.search('mar', a):
                mar_plural = re.search('\+(\d)',a)
                if mar_plural:
                    mar += int(mar_plural.groups()[0])
                    if verbose:
                        verbose_dict[a] = f'markers_added: {int(mar_plural.groups()[0])}'
                else:
                    mar += 1
                    if verbose:
                        verbose_dict[a] = 'markers_added: 1'
            
            #detecting presence on monosomies
            if re.fullmatch('-\d+|-[XY]', a):
                mono += 1
                if verbose:
                    verbose_dict[a].append('monosomy count + 1')
                #detecting number of structular changes
                 #if re.search('[inv|t].*\).*\)', a): # old definition


            #detecting presence of polysomies
            if re.fullmatch('\+\d+|\+[XY]', a):
                poly += 1
                if verbose:
                    verbose_dict[a].append('polysomy count + 1')

            #searching for structural changes
            if not re.fullmatch('[+\-]\d+c?|[+\-][XY]c?', a):
                struc += 1
                if verbose:
                    verbose_dict[a].append('Structural count + 1')
            
            if re.search('der', a):
                der += 1
                if verbose:
                    verbose_dict[a].append('der count + 1')
            
            #working out if any abnormalities are to do with 17p
            if re.search('17.*p|-17',a):
                m = re.search('\d\d?;\d\d?', a)
                if m:
                    split = re.split(';', m.group())
                    if re.findall('[pq]', a)[split.index('17')] == 'p':
                        seventeen_p = True
                        if verbose:
                            verbose_dict[a].append('17p')
                    else:
                        seventeen_p = False
                else:
                    seventeen_p = True
                    if verbose:
                        verbose_dict[a].append('17p')
            
            #looping through pre-set properties
            for p in prop_dict:
                if re.search(p, a):
                    col_true.add(p)
                    if verbose:
                        verbose_dict[a].append(prop_dict[p])
                    
    #abnormality count is equal to all non-removed + der is double counted + mar count
    abnorms = abnorms.difference(removed)
    row['Number of cytogenetic abnormalities'] = len(abnorms) + der + mar #no longer having a mar count
    
    row['Monosomy'] = mono
    row['Polysomy'] = poly
    row['Structural'] = struc
    row['abnormal(17p)'] = seventeen_p
    for c in col_true:
        row[prop_dict[c]] = True
    if verbose:
        row['verbose'] = verbose_dict
    return row


def update_using_fish(cyto_results, fish_results):
    '''
    Function that updates cyto_results using fish_results,
    overriding if discrepancies exist.
    
    Params:
    -------
    cyto_results: dict
        results from reading cytogenetic text reports
        
    fish_results: dict
        manually inserted fish results, with each result inserted
        represented by a key and its corresponding value is a boolean
        indication of the results
        
    Returns:
    --------
    cyto_results: dict
        updated cytogenetic results
    '''
    
    #dictionary mapping FISH variables to list of related cytogenetic booleans
    fish_dict = {'FISH_RUNX1-RUNX1T1': ['t(8;21)'],
 'FISH_CBFB-MYH11': ['inv(16)', 't(16;16)'],
 'FISH_PML-RARA': ['t(15;17)'],
 'FISH_Monosomy7': ['Monosomy7'],
 'FISH_Monosomy5': ['Monosomy5'],
 'FISH_del5q': ['del5q'],
 'FISH_MLL': ['t(v;11)'],
 'FISH_MECOM': ['inv(3)','t(3;3)']}
    
    #list of structural FISH abnormalities
    struc_list = ['FISH_CBFB-MYH11', 'FISH_RUNX1-RUNX1T1', 'FISH_PML-RARA', 
                  'FISH_MLL', 'FISH_MECOM']
    mono_list = ['FISH_Monosomy5','FISH_Monosomy7']
    
    #list of cytogenetic result abnormalities relating to v;11
    eleven_list = ['t(11;16)(q23.3;p13.3)', 't(2;11)', 't(9;11)', 't(6;11)', 't(10;11)']
    
    
    updated_results = cyto_results.copy()
    #assuming there is no discrepancy until found
    updated_results['FISH discrepancy'] = False
    
    
    for k in fish_results:
        if k not in fish_dict:
            #this should at least trigger a warning in the response but good enough for now
            continue
        if all([updated_results[i] != fish_results[k] for i in fish_dict[k]]):
            #updating discrepancy status
            updated_results['FISH discrepancy'] = True
            

            #updating counts if FISH true, cyto false
            if fish_results[k]:
                #updating cytogenetic result
                updated_results[fish_dict[k][0]] = True
                
                updated_results['Number of cytogenetic abnormalities'] += 1
                if k in struc_list:
                    updated_results['Structural'] += 1
                if k in mono_list:
                    updated_results['Monosomy'] += 1
            
            #updating all items that would trigger cyto result to be true,
            #to be false in accordance with FISH
            else:
                for i in fish_dict[k]:
                    updated_results[i] = False
                    
            #edge case: setting all (v;11) to False if MLL is False
            if k == 'FISH_MLL':
                if not fish_results[k]:
                    for e in eleven_list:
                        updated_results[e] = False
                
    return updated_results


def setup(abnormalities):
    """
    convert a list of abnormalities into a config for the extractor
    """
    prop_dict = properties_dict(karyotypes=None, properties=abnormalities)
    return prop_dict

def only_positive_results(result_dict):
    """
    when printing out results, only return abnormalities detected abnormalities
    """
    to_remove = []
    for k,v in result_dict.items():
        if v == 'False':
            to_remove.append(k)
    for k in to_remove:
        del result_dict[k]

    return result_dict


def extract_from_string(karyotype, prop_dict, bool_mode = 'string', fish = None,
                        verbose=False):
    """
    Run extraction on a single karyotype string, extraction based on prop_dict
    prop_dict can be created with setup()
    This function guarantees the output will have a property key for every abnormality value in prop_dict
    as well as some additional created by parse_karyotype
    Anything in prop_dict that is not detected will default to False
    bool_mode: return type for boolean values. If 'string' then True -> "True" and Talse -> "False", otherwise return bool. 
    """
    input = {
        'Cytogenetics': karyotype.strip(),
        'Error': False,
        'Error description': []
    }
    result = parse_karyotype(input, prop_dict, verbose=verbose)
    for abn in prop_dict.values():
        if abn not in result:
            result[abn] = False
            
    #using FISH
    if fish:
        result = update_using_fish(result, fish)
    else:
        result['FISH discrepancy'] = False

    if bool_mode == 'string':
        for abn in result:
            if abn == 'Error':
                continue
            if type(result[abn]) == bool:
                result[abn] = str(result[abn])
    output = {'error': result['Error'],
              'error_message': result['Error description'],
              'Warnings': result['Warnings'],
              'result': only_positive_results(result),
              'fish_available':False}
    
    if fish:
        output['fish_available'] = True

    if verbose:
        output['verbose'] = result['verbose']
        del result['verbose']
    
    return output

def base_extraction():
    ex = ["-Y", "-X", 'del11q', 'del12p', 
    'del13q', 'del5q', 'del7q', 'idic(X)(q13)', 'isochromosome17q', 'Monosomy13', 
    'Monosomy17', 'Monosomy5', 'Monosomy7', 't(1;3)', 't(11;16)(q23.3;p13.3)', 
    't(12p)', 't(17p)', 't(2;11)', 't(3;21)', 't(3;5)', 't(5;10)', 't(5;12)', 
    't(5;17)', 't(5;7)', 't(5q)', 't(1;22)', 'inv(3)', 't(3;3)', 't(6;9)', 
    't(9;22)', 't(16;16)', 'inv(16)', 't(8;21)', 't(15;17)', 't(9;11)', 
    't(6;11)', 't(10;11)','t(v;11)']
    return ex

def extraction_2022():
    # additionally checks for 't(8;16)(p11;p13)', 't(3q26.2;v)'
    ex = ["-Y", "-X", 'del11q', 'del12p', 
    'del13q', 'del5q', 'del7q', 'idic(X)(q13)', 'isochromosome17q', 'Monosomy13', 
    'Monosomy17', 'Monosomy5', 'Monosomy7', 't(1;3)', 't(11;16)(q23.3;p13.3)', 
    't(12p)', 't(17p)', 't(2;11)', 't(3;21)', 't(3;5)', 't(5;10)', 't(5;12)', 
    't(5;17)', 't(5;7)', 't(5q)', 't(1;22)', 'inv(3)', 't(3;3)', 't(6;9)', 
    't(9;22)', 't(16;16)', 'inv(16)', 't(8;21)', 't(15;17)', 't(9;11)', 
    't(6;11)', 't(10;11)', 't(8;16)(p11;p13)', 't(3q26.2;v)',
          't(v;11)']
    return ex

def available_configs():
    configs = {}
    configs["BJH2021"] = setup(base_extraction())
    configs["ELN2022"] = setup(extraction_2022())
    return configs

VERBOSE = True
if __name__ == '__main__':
    # process from excel file
#    karyotypes = load_file()
#    prop_dict = properties_dict(karyotypes=karyotypes)
#    karyotypes['Cytogenetics'] = karyotypes.apply(remove_artefact, axis=1)
#    karyotypes['Error'] = False
#    karyotypes['Error description'] = None
#    results = karyotypes.apply(parse_karyotype, axis=1, args=(prop_dict,))
#    results.loc[results['Error']==False] = results.loc[results['Error']==False].fillna(False)
#    results['Error'] = results['Error'].astype(bool)
#    results.to_excel('Cytogenetics_output_V4.xlsx')

    #process from single string as in API
    #abn = base_extraction()
    abn = extraction_2022()
    props = setup(abn)
    fish_results = {'FISH_RUNX1-RUNX1T1': False,
     'FISH_CBFB-MYH11': False,
     'FISH_PML-RARA': False,
     'FISH_Monosomy7': False,
     'FISH_Monosomy5': False,
     'FISH_del5q': False,
     'FISH_MLL': False,
     'FISH_MECOM': False
    }
    #report = "  45,X,-Y[17]/46,XY[3]   "
    #report = "47,XY,+21c[6]/48,idem,+11,der(19)t(1;19)(q23;p13.3)[4]"
    #report = "47,XY,+11,t(3;19)(q26.2;p13.3)[4]"
    #report = " 46,xx,t(8;16)(p11.2;p13.3)[20]"
    #report = "46,XY, -17[16]"
    #report = "45,XX,t(3;21)(q26;q?11.2),del(5)(q23-31q33),-7[14]"
    report = "46,XX,t(3;12)(q26;p13)[18]"
    result = extract_from_string(report, props, fish=fish_results, verbose = VERBOSE)
    print(report)
    print(result)
    #print(props)