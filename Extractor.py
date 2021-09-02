#importing modules

import re
import pandas as pd
#import streamlit as st
#import base64

punct = r'()[]/'
punct_dict = {s:0 for s in punct}


def load_file(file= 'Cytogenetics_TS_Apr2021.xlsx', streamlit=False):
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
    if streamlit:
        file = st.file_uploader("Upload excel here")
    karyotypes = pd.read_excel(file)    
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
            f = v[:m.start()] + '(' + m.group()[:-1] + \
                ')(' + m.group()[-1]+ v[m.end():]
            v=f
            
        #creating escape characters for strings
        v = re.escape(v)
        
        #creating special case for for t(v;11)
        if v == 't\\(v;11\\)':
            v = '(t\\(\d+;11\\))|(t\\(11;\d+\\))'
        
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


def gram_error(string):
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
    
    #basic report errors
    if string == 'Error':
        error.append('String report missing')
        
    if re.search('fail', string.lower()):
        error.append('String report indicates failure')
    
    #highlight if error contains uncertainty
    if re.search('\\?', string):
        error.append('Question mark in string means report is uncertain')
        
    #all reports should include commas somewhere 
    if not re.search(',', string):
        error.append('Missing comma')
        
    #checking for constitutional changes
    if re.search('[^a-z]?c[^p]', string) or re.search('[^a-z]c$', string):
        error.append('constitutional changes present')
    
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
                    error.append(f'Variable number of markers in report. Minimum number ({d}) used by default')
        
        #in the case of a range of potential chromosome counts
        if re.search('[~]', chrom):
            try:
                low_num, high_num = int(chrom[:2]), int(chrom[-2:])
            except ValueError:
                error.append('Part of report not clearly defined by two chromosome numbers followed by comma (e.g. "43~45,")')
                continue
            if low_num <= expected <= high_num:
                pass
            elif expected > high_num:
                error.append(f'chromosome number lower than expected in subsection number {i+1}')
            else:
                error.append(f'chromosome number higher than expected in subsection number {i+1}')
                
        #if only one chromosome count is present
        else:
            try:
                num = int(chrom[:2])
            except ValueError:
                error.append('Start of report missing clear chromosome number followed by comma (e.g. "46,")')
                continue
            if expected == num:
                pass
            elif expected < num:
                error.append(f'chromosome number higher than expected in subsection number {i+1}')
            else:
                error.append(f'chromosome number lower than expected in subsection number {i+1}')
    return error

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

def check_trans_dict(trans_dict, col_true):
    '''
    Takes in multi translocation dictionary and detects abnormalities associated
    with adjacent chromosome translocation
    '''
    chr_keys = list(trans_dict)
    for i in range(len(chr_keys)-1):
        sorted_list = sorted([chr_keys[i],chr_keys[i+1]])
        first_chr, second_chr = sorted_list[0], sorted_list[1]
        for p in prop_dict:
            if re.search(p, f't({first_chr};{second_chr})'):
                col_true.add(p)
    return col_true

def parse_karyotype(row, prop_dict):
    '''
    Working row by row on a dataframe, detects abnormalities
    and fills in boolean and integer counts.
    '''
    error = gram_error(row['Cytogenetics'])
    permissible_errors = ['chromosome', 'Variable', 'constitutional', 'Question']
    
    #checking for error
    if error:
        row['Error description'] = error
        for e in error:
            for p_e in permissible_errors:
                if e.startswith(p_e):
                    break
            else:
                row['Error'] = True
                return row
    
    #initialising counts and sets for abnormalities
    abnorms = set(re.split('/|,|\[', row['Cytogenetics']))
    removed = set()
    col_true = set()
    mono = 0
    struc = 0
    der = 0
    mar = 0
    seventeen_p = False
    for a in abnorms:
        #cremoving normal report segments
        if re.fullmatch('(\d\d([~-]\d\d)?|X[XY]?|(cp)?\d\d?\]|idem)(\??c)?', a):
            removed.add(a)
#         if re.search('mar', a): #no longer the case to remove mar
#             removed.add(a)

        #checking abnormalities
        if a not in removed:
            
            #checking multi-translocations
            if re.search('t\(\d\d?;\d\d?;\d\d?',a):
                trans_dict = make_multi_translocation_dict(a)
                col_true = check_trans_dict(trans_dict, col_true)
            
            #counting number of markers present
            if re.search('mar', a):
                mar_plural = re.search('\+(\d)',a)
                if mar_plural:
                    mar += int(mar_plural.groups()[0])
                else:
                    mar += 1
            
            #detecting presence on monosomies
            if re.fullmatch('-\d+|-[XY]', a):
                mono += 1
            #detecting number of structular changes
            if re.search('[inv|t].*\).*\)', a):
                struc += 1
            
            if re.search('der', a):
                der += 1
            
            #working out if any abnormalities are to do with 17p
            if re.search('17.*p|-17',a):
                m = re.search('\d\d?;\d\d?', a)
                if m:
                    split = re.split(';', m.group())
                    if re.findall('[pq]', a)[split.index('17')] == 'p':
                        seventeen_p = True
                    else:
                        seventeen_p = False
                else:
                    seventeen_p = True
            
            #looping through pre-set properties
            for p in prop_dict:
                if re.search(p, a):
                    col_true.add(p)
                    
    #abnormality count is equal to all non-removed + der is double counted + mar count
    abnorms = abnorms.difference(removed)
    row['Number of cytogenetic abnormalities'] = len(abnorms) + der + mar
    
    row['Monosomy'] = mono
    row['Structural'] = struc
    row['abnormal(17p)'] = seventeen_p
    for c in col_true:
        row[prop_dict[c]] = True
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
        if all([updated_results['result'][i] != fish_results[k] for i in fish_dict[k]]):
            #updating discrepancy status
            updated_results['FISH discrepancy'] = True
            

            #updating counts if FISH true, cyto false
            if fish_results[k]:
                #updating cytogenetic result
                updated_results['result'][fish_dict[k][0]] = True
                
                updated_results['result']['Number of cytogenetic abnormalities'] += 1
                if k in struc_list:
                    updated_results['result']['Structural'] += 1
                if k in mono_list:
                    updated_results['result']['Monosomy'] += 1
            
            #updating all items that would trigger cyto result to be true,
            #to be false in accordance with FISH
            else:
                for i in fish_dict[k]:
                    updated_results['result'][fish_dict[k][i]] = False
                    
            #edge case: setting all (v;11) to False if MLL is False
            if k == 'FISH_MLL':
                if not fish_results[k]:
                    for e in eleven_list:
                        updated_results['result'][e] = False
                
    return updated_results


def setup(abnormalities):
    """
    convert a list of abnormalities into a config for the extractor
    """
    prop_dict = properties_dict(karyotypes=None, properties=abnormalities)
    return prop_dict

def extract_from_string(karyotype, prop_dict, bool_mode = 'string', fish = None):
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
        'Error description': ""
    }
    result = parse_karyotype(input, prop_dict)
    for abn in prop_dict.values():
        if abn not in result:
            result[abn] = False
    if bool_mode == 'string':
        for abn in result:
            if abn == 'Error':
                continue
            if type(result[abn]) == bool:
                result[abn] = str(result[abn])
    output = {'error': result['Error'], 'error_message': result['Error description'], 'result': result, 'fish_available':False}
    
    #using FISH
    if fish:
        output['fish_available'] = True
        output = update_using_fish(output, fish)
    return output

def base_extraction():
    ex = ["-Y", "-X", 'del11q', 'del12p', 
    'del13q', 'del5q', 'del7q', 'idic(X)(q13)', 'isochromosome17q', 'Monosomy13', 
    'Monosomy17', 'Monosomy5', 'Monosomy7', 't(1;3)', 't(11;16)(q23.3;p13.3)', 
    't(12p)', 't(17p)', 't(2;11)', 't(3;21)', 't(3;5)', 't(5;10)', 't(5;12)', 
    't(5;17)', 't(5;7)', 't(5q)', 't(1;22)', 'inv(3)', 't(3;3)', 't(6;9)', 
    't(9;22)', 't(16;16)', 'inv(16)', 't(8;21)', 't(15;17)', 't(9;11)', 
    't(6;11)', 't(10;11)', 't(v;11)']
    return ex


if __name__ == '__main__':
    # process from excel file
    karyotypes = load_file()
    prop_dict = properties_dict(karyotypes=karyotypes)
    karyotypes['Cytogenetics'] = karyotypes.apply(remove_artefact, axis=1)
    karyotypes['Error'] = False
    karyotypes['Error description'] = None
    results = karyotypes.apply(parse_karyotype, axis=1, args=(prop_dict,))
    results.loc[results['Error']==False] = results.loc[results['Error']==False].fillna(False)
    results['Error'] = results['Error'].astype(bool)
    results.to_excel('Cytogenetics_output_V4.xlsx')

    #process from single string as in API
    abn = base_extraction()
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
    report = "  45,X,-Y[17]/46,XY[3]   "
    result = extract_from_string(report, props, fish=fish_results)
    print(report)
    print(result)