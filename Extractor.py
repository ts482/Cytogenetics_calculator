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
    #karyotypes = karyotypes.drop(columns='ID')
    
    return karyotypes

def dotreplace(matchobj):
    '''
    helper function for re.sub
    replaces chromosome decimals (e.g. "17p14(.2)"), with a string that makes
    these decimal points optional
    
    LEGACY COMMENTED OUT CODE - accepts +/- .1 of the chromosome
    '''
    #int_match = int(matchobj.group(2))
    #first_num = int_match-1
    #second_num = int_match+1
    #first_num, second_num = str(first_num), str(second_num)
    #return f'{matchobj.group(1)}(?:\.[{first_num}-{second_num}])?'
    return f'{matchobj.group(1)}(?:\.{matchobj.group(2)})?'



def substitute_v(t):
    #important to note f curly braces and regex curly braces
    der = re.search('i?der\(([0-9XxYy]{1,2})\)t', t)
    if der:
        der_chrom = der.group(1)
    m = re.search('v;([0-9XxYy]{1,2})([qp]\d{0,2}(?:\.\d)?)?', t)
    if not m:
        m = re.search('([0-9XxYy]{1,2})([qp]\d{0,2}(?:\.\d)?)?;v', t)
    chrom = m.group(1)
    arm = m.group(2)
    if arm:
        armdot = re.search('\.(\d)', arm)
        armnum = re.search('[qp]\d{1,2}', arm)
        if armdot:
            dotnum = armdot.group(1)
            #first_num = str(int(dotnum)-1) #legacy code to do XpY.Z, where Z+/-1
            #second_num = str(int(dotnum)+1)
            #arm = re.sub('\.\d+', f'(?:\.[{first_num}-{second_num}])?', arm)
            arm = re.sub('\.\d+', f'(?:\.{dotnum})?', arm)

        elif armnum:
            arm = arm + '(?:\.\d)?'
        elif not armnum:
            arm = arm + '(?:\d{1,2})?(?:\.\d)?'
        if der:
            return f'(i?der\({der_chrom}\)' + 't\((?:[0-9XxYy]){1,2};' + \
                f'{chrom}\)' + '\([qp]\d{1,2}(?:\.\d)?;' + f'{arm}\))' + \
                f'|(i?der\({der_chrom}\)t\({chrom};' + '(?:[0-9XxYy]){1,2}\)' + \
                f'\({arm};' + '[qp]\d{1,2}(?:\.\d)?\))'
        return 't\((?:[0-9XxYy]{1,2});' + f'{chrom}\)' + \
            '\([qp]\d{1,2}(?:\.\d)?;' + f'{arm}\)' + \
            f'|t\({chrom};' + '(?:[0-9XxYy]){1,2}\)'+ f'\({arm};' + \
                '[qp]\d{1,2}(?:\.\d)?\)'
    if der:
        return f'(i?der\({der_chrom}\)' + 't\((?:[0-9XxYy]){1,2};' + \
            f'{chrom}\))|(i?der\({der_chrom}\)t\({chrom};' + '(?:[0-9XxYy]){1,2}\))'
    return '(t\((?:[0-9XxYy]){1,2};' + f'{chrom}\))|(t\({chrom};' + \
        '(?:[0-9XxYy]){1,2}\))'

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
        
        #formatting trisomies
        if v.startswith('Trisomy'):
            v = '+' + v[7:]
        
        #replacing all chromosome dots with contingent verisions which are +/- .1
        
        
        if not re.search('\(v;|;v\)', v):
            if not re.search('inv',v): #this bracket insertion is not applicable to inversions
                #correctly formatting translocations by seperating p and q
                matches = re.finditer('([0-9XxYy]{1,2})(p|q)',v)
                for m in matches:
                    #different rule for translocations - removing first bracket
                    if re.search('(?:i|t|inv|del|add)\(', m.string):
                        f = v[:m.start()] + m.group()[:-1] + \
                            ')(' + m.group()[-1]+ v[m.end():]
                    else:
                        f = v[:m.start()] + '(' + m.group()[:-1] + \
                            ')(' + m.group()[-1]+ v[m.end():]        
                    v=f
        
        
        
            #creating escape characters for strings
            v = re.escape(v)
            
            #formatting exception for deletions and additions
            if re.search('^(add|del|i)', v):
                v = re.sub('\\\\\)\\\\\(', '(?:\)\()?', v)
                v = re.sub('\\\\\)$', '', v)
            
            #replacing dots with less specific versions
            if re.search('[pq]\d{1,2}\\\.\d',v):
                v = re.sub('([pq]\d{1,2})\\\.(\d)', dotreplace, v)
                
            #putting brackets at the end to stop catching
            #wrong number on second chromosomes
            if re.search('^(del|add|idic|i|inv)',v):
                porq = re.search('([pq])$',v)
                if porq:
                    v = v + '\d{1,2}?(?:\.\d)?(?:[pq]\d{1,2}?(?:\.\d)?)?|' + v[:-1] + \
                        f'[pq]\d?\d?(?:\.\d)?{porq.group(1)}\d?\d?(?:\.\d)?'
                else:
                    v = v + '\)'
                
        else:
            v = substitute_v(v)
        
        #special case for DEK/NUP214 fusion
        if v == 't\(6;9\)\(p22(?:\.3)?;q34(?:\.1)?\)':
            v = 't\(6;9\)\(p2(?:3|2(?:\.3)?);q34(?:\.1)?\)'
        
        #making sure monosomies and trisomies are counted from start:
        if re.search('\\\[-+][1-9XxYy]{1,2}', v):
            v = re.sub('\\\\', '^\\\\', v)
        
        #special case for dots
        #if v == "t\(8;16\)\(p11;p13\)":
        #    v = "t\(8;16\)\(p11\.?\d?;p13\.?\d?\)"
        
        #creating special case for for t(v;11)
        #if v == 't\\(v;11\\)':
        #    v = '(t\\((?:[0-9XxYy]+);11\\))|(t\\(11;(?:[0-9XxYy]+)\\))'
            
        #another special case for t(3q26.2;v)
        #if v == 't\\(3\\)\\(q26\\.2;v\\)':
        #    v = '(t\((?:[0-9XxYy]+);3\)\([pq]\d+(\.\d+)?;q26(\.2)?\))|(t\(3;(?:[0-9XxYy]+)\)\(q26(\.2)?;[pq]\d+(\.\d+)?\))'

        
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
    
    
    #ERRORS
    
    #basic report errors
    if string == 'Error':
        error.append('String report missing')
        
    if re.search('fail', string.lower()):
        error.append('String report indicates failure')
        
    if re.search('or', string.lower()):
        error.append('String contains \'or\' and therefore cannot be interpreted')
    
    if re.search('^\D', string):
        error.append('Report must start with a number')
        
    if not re.search('\]$', string):
        error.append('report must end in \']\' ')

    #all reports should include commas somewhere 
    if not re.search(',', string):
        error.append('Missing comma')

    #incorrectly counts abnormalities if split character are together
    if re.search('[/,[]{2}', string):
        if re.search('//', string):
            error.append('The karyotype appears to be a chimeric karyotype post stem cell transplant. The app is not designed to analyse these reports')
        else:
            error.append('The following symbols must not be next to another or themselves: / , [')
        
    if re.search('(?:p|q)ter|::', string):
        seq = re.search('(?:p|q)ter|::', string).group(0)
        error.append(f"'{seq}' - this extractor has not been set up to detect karyotypes written in this format")
    
    ## WARNINGS
    
    #multiple clones contain markers
    if re.search('mar.*/.*mar', string):
        warning.append('multiple clones contain markers, marker counts may not be accurate')
    
    #checking for YX

    if re.search('YX', string):
        warning.appened('"YX" appears to be wrong way round. Did you mean "XY"?')
        
    #checking for constitutional changes
    if re.search('[^a-z]c[^p]', string) or re.search('[^a-z]c$', string):
        warning.append('constitutional changes present')
        
    if re.search('<\dn>', string):
        warning.append('polyploidy present')
        
    #highlight if error contains uncertainty
    if re.search('\\?', string):
        warning.append('Question mark in string means report is uncertain')
        
    if re.search('inc', string):
        warning.append('Note incomplete karyotype reported, additional abnormalities may not have been identified by this test')
        
    if re.search('t\([0-9XxYy]{1,2};[0-9XxYy]\d?;[0-9XxYy]\d?',string):
        warning.append('translocation with more than 2 parts detected. Results should be treated with caution.')
    
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
        if not re.search('(idem)|(sd?l)', s):
            expected = 46
        else:
            if re.search('(idem)|(sd?l\d?(x\d)?)', s):
                warning.append('stemline duplication present, results should be taken with caution')
        expected -= len(re.findall('\-', s))
        expected += len(re.findall('\+', s))
        if re.search('mar', s):
            mar_plural = re.search('\+(\d~)(\d)?mar',s)
            if mar_plural:
                expected += int(mar_plural.groups()[1]) - 1
                tilda_present =  re.search('\+(\d)~(\d)mar', s)
                #error if variable number of markers included (e.g. +2~4mar)
                if tilda_present:
                    d = tilda_present.group(2)
                    warning.append(f'Variable number of markers in report. Maximum number ({d}) used by default')
        if verbose:
            chr_count['Expected number'] = expected
        
        #in the case of a range of potential chromosome counts
        if re.search('[~-]', chrom):
            try:
                low_num, high_num = map(lambda x: int(x),re.search('(\d+)[~-](\d+)',chrom).groups([1,2]))
            except ValueError:
                error.append('Part of report not clearly defined by two chromosome numbers followed by comma (e.g. "43~45,")')
                continue
            
            warning.append(f'variable number of chromosomes detected in subsection {i+1}. Higher number used.')
            
            if high_num > 64:
                warning.append(f'high chromosome number detected indicating polyploidy: {high_num}. Abnormality counts should be treated with caution.')
                while expected + 18 < high_num:
                    expected += 23
                
            elif high_num > 49:
                warning.append(f' high chromosome number detected indicating hyperploidy: {high_num}. Abnormality counts should be treated with caution.')
        
            if verbose:
                chr_count['Reported number'] = str(low_num) + '~' + str(high_num)
            if low_num <= expected <= high_num:
                pass
            elif expected > high_num:
                warning.append(f'chromosome number lower than expected in subsection number {i+1}')
            else:
                warning.append(f'chromosome number higher than expected in subsection number {i+1}')
                
        #if only one chromosome count is present
        elif not error: #making sure no crashes occur from earlier spotted error
            try:
                num = int(re.search('^(\d+)', chrom).group(0))
            except ValueError:
                error.append('Start of report missing clear chromosome number followed by comma (e.g. "46,")')
                continue
            if num > 64:
                warning.append(f'high chromosome number detected indicating polyploidy: {num}.')
                while expected + 18 < num:
                    expected += 23
            elif num > 49:
                warning.append(f' high chromosome number detected indicating hyperploidy: {num}')
            if verbose:
                chr_count['Reported number'] = num
            if expected == num:
                pass
            elif expected < num:
                warning.append(f'chromosome number higher than expected in subsection number {i+1}')
            else:
                warning.append(f'chromosome number lower than expected in subsection number {i+1}')
        #checking for composite clones
        if len(warning)>0:
            if re.search('chromosome', warning[-1]):
                if re.search('cp', string):
                    cp_cr = warning.pop()
                    cp_cr = cp_cr + ', but note this is a composite clone'
                    warning.append(cp_cr)
    return error, warning, chr_count

def base_translocation_check(translocation, col_true, prop_dict,
                             verbose=False):
    '''
    

    Parameters
    ----------
    translocation : string
        
    col_true : dict
        DESCRIPTION.
    prop_dict : dict
        DESCRIPTION.
    verbose : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    '''
    
    
    
    
def make_multi_translocation_dict(string):
    '''
    Creates a translocation list if only full chromosomes are involved
    in translocation, or a dictionary if chromosome arms are involved
    '''
    
    #finding out how many chromosomes in translocation
    exp = '([0-9XxYy]{1,2});([0-9XxYy]{1,2})'
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

def check_trans_dict(trans_dict, col_true,prop_dict, verbose=False):
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
        if isinstance(trans_dict, dict):
            first_arm, second_arm = trans_dict[first_chr], trans_dict[second_chr]
        for p in prop_dict:
            if isinstance(trans_dict, list):
                if re.search(p, f't({first_chr};{second_chr})'):
                    col_true.add(p)
                    if verbose:
                        verbose_list.append(p)
            elif isinstance(trans_dict, dict):
                if re.search(p, f't({first_chr};{second_chr})({first_arm};{second_arm})'):
                    col_true.add(p)
                    if verbose:
                        verbose_list.append(prop_dict[p])
    if verbose:
        return col_true, verbose_list
    return col_true

def two_part_translocation_check(two_part_translocation,col_true, prop_dict,
                                 verbose=False):
    if verbose:
        verbose_list = []
    match = re.search('t\(([XxYy\d]{1,2});([XxYy\d]{1,2})\)\(([pq])\d*(?:\.\d+)?;([pq])\d*(?:\.\d+)?\)',
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

def split_loci_into_num_and_arm(locus):
    '''
    Helper function for setup warning dict

    Parameters
    ----------
    loci : STR
        a string representing a loci (that is likely clinically relevant)

    Returns
    -------
    dictionary that splits the loci number from the chromsome arm

    '''
    
    num, arm = re.search('(\d+)([qp].*)',locus).group(1,2)
    return {'num':num, 'arm':arm}

def setup_suspicious_dict():
    loci = ['11q23.3', '3q26.2','11p15.4','17q21.2']
    sus_dict = {l:split_loci_into_num_and_arm(l) for l in loci}
    return sus_dict

def parse_karyotype_clone(row, prop_dict, verbose=False):
    '''
    Working row by row on a dataframe, detects abnormalities
    and fills in boolean and integer counts.
    Second function version: using information about each clone
    '''
    error, warning, chr_count = gram_error(row['Cytogenetics'])
    if verbose:
        row['chr_count'] = chr_count
        
    
    #checking for error
    if len(error)>0:
        row['Error description'] = error
        row['Error'] = True
        if verbose:
            row['verbose'] = "see error description"
        return row
    
    

    #preparing string before splitting to abnormalities
    full_string = row['Cytogenetics']
    full_string = re.sub('\?', '', full_string)
    
    #splitting into clones
    clones = list(re.split('/', full_string))
    
    #clone_dict = {}
    sl_dict = {}
    
    all_clone_abnorms = []
    abnorm_set = set()
    
    removed_string = '(\d\d([~-]\d\d)?(<\dn>)?|[XxYy]+|inc|(cp)?\d\d?\]|idem|sd?l\d?)(\??c)?(x\d+)?'
    
    for i,clone in enumerate(clones):
    
        #initialising counts and sets for abnormalities
        clone_abnorms = list(re.split(',|\[', clone))
        
        
        
        plicates = list()
        
        removed = set()
        
        idem_sl_multiplier = 1
        for a in clone_abnorms:
            multiplier = re.search('x(\d)$', a)
            if multiplier:
                if re.search('idem|sd?l\d?',a):
                    idem_sl_multiplier = int(multiplier.group(1))
                else:
                    non_x_string = re.sub('x(\d)$', '', a)
                    for i in range(int(multiplier.group(1))): 
                        plicates.append(non_x_string)
                    removed.add(a)
                    
        clone_abnorms = clone_abnorms + plicates
                    
                
        if re.search('idem',clone):
            clone_abnorms = clone_abnorms + all_clone_abnorms[0] * idem_sl_multiplier
        
        #looking for presence of sl groups
        sl_group = re.search('(sd?l\d?)', clone)
        if sl_group:
            if sl_group.group(0) in sl_dict:
                clone_abnorms = clone_abnorms + sl_dict[sl_group.group(0)] * idem_sl_multiplier
            else:
                clone_abnorms = clone_abnorms + all_clone_abnorms[i-1] * idem_sl_multiplier
                sl_dict[sl_group.group(0)] = all_clone_abnorms[i-1]
        
       
        
        for a in clone_abnorms:
            if re.fullmatch(removed_string,a):
                removed.add(a)
            #removing constitutional changes
            if re.search('[XxYy1-9]c$', a):
                removed.add(a)
        
        clone_abnorms = [ab for ab in clone_abnorms if ab not in removed]
        
        #finding out abnormalities which cancel out
        for i, ab in enumerate(clone_abnorms):
            somy = re.search('(^[+-])(.*)', ab)
            if somy:
                if somy.group(1) == '+':
                    for ab2 in clone_abnorms:
                        if ab2 == '-' + somy.group(2):
                            clone_abnorms.remove(ab)
                            clone_abnorms.remove(ab2)
                            break
                elif somy.group(1) == '-':
                    for ab2 in clone_abnorms:
                        if ab2 in [somy.group(2),'+' + somy.group(2)]:
                            clone_abnorms.remove(ab)
                            clone_abnorms.remove(ab2)
                            break

        for ab in clone_abnorms:
            abnorm_set.add(ab)
        all_clone_abnorms.append(clone_abnorms)
    
    
    abnorm_list = list(abnorm_set)
    abnorm_count = pd.DataFrame(index = clones, columns = abnorm_list)
    for clone, abnorm_l in zip(clones, all_clone_abnorms):
        for unique_abn in abnorm_set:
            abnorm_count.loc[clone, unique_abn] = abnorm_l.count(unique_abn)
        
    
    final_abn_count = abnorm_count.max(axis=0)
    abnorms = []
    for abn, count in final_abn_count.items():
        abnorms = abnorms + [abn] * count
    
    
    
    verbose_dict = {}
    
    col_true = set()
    mono = 0
    non_sex_mono = 0
    poly = 0
    struc = 0
    seventeen_p = False
    der = 0
    mar = 0
    er_mar = 0 #counts the erroneously cytogenetic count for markers
    for a in abnorms:
        if verbose:
            verbose_dict[a] = [f'abnormality count = {final_abn_count[a]}']
        #removing normal report segments
        if re.fullmatch(removed_string, a):
            removed.add(a)
            if verbose:
                if a.endswith(']'):
                    del verbose_dict[a]
                else:
                    verbose_dict[a] = 'Normal'


        #checking abnormalities
        if a not in removed:
            
            #checking multi-translocations
            if re.search('t\([0-9XxYy]{1,2};[0-9XxYy]\d?;[0-9XxYy]\d?',a):
                trans_dict = make_multi_translocation_dict(a)
                if verbose:
                    col_true, verbose_list = \
                        check_trans_dict(trans_dict, col_true,prop_dict,
                                         verbose=True)
                    for v in verbose_list:
                        verbose_dict[a].append(v)
                else:
                    col_true = check_trans_dict(trans_dict, prop_dict, col_true)

            #checking two-part translocations
            if re.search('t\(([XxYy\d]{1,2});([XxYy\d]{1,2})\)\(([pq])\d*(\.\d+)?;([pq])\d*(\.\d+)?\)',
                         a):
                if verbose:
                    col_true, verbose_list = two_part_translocation_check(a,
                        col_true, prop_dict, verbose=True)
                    for v in verbose_list:
                        verbose_dict[a].append(v)
                else:
                    col_true = two_part_translocation_check(a,
                                col_true, prop_dict, verbose-False)
                
            
            #detecting presence on monosomies
            if re.fullmatch('-[0-9]{1,2}', a): #used to be [0-9XxYy]
                mono += 1
                if re.fullmatch('-[0-9]{1,2}', a):
                    non_sex_mono +=1
                    if verbose:
                        verbose_dict[a].append('non sex chromosome monosomy count + 1')
                if verbose:
                    verbose_dict[a].append('monosomy count + 1')
                #detecting number of structular changes
                 #if re.search('[inv|t].*\).*\)', a): # old definition


            #detecting presence of polysomies
            if re.fullmatch('\+[0-9XxYy]{1,2}', a):
                poly += 1
                if verbose:
                    verbose_dict[a].append('polysomy count + 1')

            #searching for structural changes
            if not re.fullmatch('[+\-][0-9XxYy]{1,2}c?', a):
                #markers are also structural but counted elsewhere, higher in the code
                if not re.search('mar',a):
                    struc += 1
                    if verbose:
                        verbose_dict[a].append('Structural count + 1')
            
            if re.search('i?der', a):
                der += 1
                if verbose:
                    verbose_dict[a].append('der count + 1')
                    
                    
            #counting number of markers present
            if re.search('mar', a):
                er_mar += 1
                mar_plural = re.search('\+(\d)',a)
                if mar_plural:
                    variable_markers = re.search('\+\d[~-](\d)', a)
                    if variable_markers:
                        mar_added = int(variable_markers.group(1))
                    else:
                        mar_added = int(mar_plural.group(1))
                    mar += mar_added
                    if verbose:
                        verbose_dict[a].append(f'markers_added: {mar_added}')
                else:
                    mar += 1
                    if verbose:
                        verbose_dict[a].append('markers_added: 1')
            
            #looping through pre-set properties
            for p in prop_dict:
                if re.search(p, a):
                    col_true.add(p)
                    if verbose:
                        verbose_dict[a].append(prop_dict[p])
            
            #working out if any abnormalities are to do with 17p
            if any([re.search('17(?:\)\()?p',a), 
                    re.search('-17',a), 
                    re.search( 't\\((?:[0-9XxYy]{1,2});17\\)\\([qp]\\d{1,2}(?:\\.\\d)?;'+
                              'p(?:\\d{1,2})?(?:\\.\\d)?\\)|'+
                              't\\(17;(?:[0-9XxYy]){1,2}\\)\\(p(?:\\d{1,2})?(?:\\.\\d)?;'+
                              '[qp]\\d{1,2}(?:\\.\\d)?\\)', a)]):
                seventeen_p = True
                if verbose:
                    verbose_dict[a].append('17p')
                    
            #adding warnings onto suspicious abnormalities:
            sus_dict = setup_suspicious_dict()
            for k,v in sus_dict.items():
                sus_num = v['num']
                sus_arm = v['arm']
                if re.search(f'{sus_num}(?:\)\()?{sus_arm}',a):
                    warning.append(f'mutation found: {k}. This is not a' + \
                        ' pre-defined mutation, but it may be clinically relevant')
                    if verbose:
                        verbose_dict[a].append(f'suspicious mutation: {k}')
            
                    
    #abnormality count is equal to all non-removed + der is double counted
    row['Number of cytogenetic abnormalities'] = len(abnorms) + der + mar - er_mar

    row['Monosomy'] = mono
    row['Polysomy'] = poly
    row['Structural'] = struc + mar + der
    row['abnormal(17p)'] = seventeen_p
    
    row['NonSexChromosomeMonosomies'] = non_sex_mono
    
    #including warnings
    row['Warnings'] = warning
    
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
                        verbose=False, only_positive=False):
    """
    Run extraction on a single karyotype string, extraction based on prop_dict
    prop_dict can be created with setup()
    This function guarantees the output will have a property key for every abnormality value in prop_dict
    as well as some additional created by parse_karyotype
    Anything in prop_dict that is not detected will default to False
    bool_mode: return type for boolean values. If 'string' then True -> "True" and Talse -> "False", otherwise return bool. 
    """
    input_dict = {
        'Cytogenetics': re.sub('\s', '',karyotype),
        'Error': False,
        'Error description': []
    }
    result = parse_karyotype_clone(input_dict, prop_dict, verbose=verbose)
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
    if only_positive:
        result = only_positive_results(result)
    
    output = {'error': result['Error'],
              'error_message': result['Error description'],
              'Warnings': result['Warnings'],
              'result': result,
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
    ex = [ # sex chromosome monosmies 
    "-Y", "-X",
    #deletions
    'del5q', 'del7q', 'del11q', 'del12p', 'del13q',
    'del(17)(q21.2q21.2)', 'del(17p)', 'del(20q)',
    #additions
    'add(5q)','add(7q)', 'add(12p)', 'add(17p)',
    #inversions
    'idic(X)(q13)', 'i(17q)', 'inv(3)(q21.3q26.2)','inv(16)(p13.1q22)',
    'inv(16)(p13.3q24.3)',
    #monosomies
    'Monosomy5', 'Monosomy7', 'Monosomy13', 'Monosomy17',
    #trisomies
    'Trisomy8',
    #translocations ordered by first chromosome, then second chromosome:
    #t1
    't(1;3)(p36.3;q21.3)','t(1;22)(p13.3;q13.1)',
    #t2
    't(2;11)(p21;q23.3)', 
    #t3
    't(3;3)(q21.3;q26.2)', 't(3;5)(q25.3;q35.1)', 't(3;21)(q26.2;q22.1)', 't(3q26.2;v)',
    #t5
    't(5;7)(q32;q11.2)', 't(5;10)(q32;q21)', 't(5;11)(q35.2;p15.4)',
    't(5;12)(q32;p13.2)', 't(5;17)(q32;p13.2)', 't(v;5q)', 
    #t6
    't(6;9)(p22.3;q34.1)', 't(6;11)',
    #t7
    't(7;12)(q36.3;p13.2)',
    #t8
    't(8;16)(p11.2;p13.3)', 't(8;21)(q22;q22.1)',
    #t9
    't(9;11)(p21.3;q23.3)', 't(9;22)(q34.1;q11.2)',
    #t10
    't(10;11)(p12.3;q14.2)',
    #t11
    't(11;12)(p15.4;p13.3)', 't(11;16)(q23.3;p13.3)',
    't(v;11p15.4)', 't(v;11q23.3)',
    #t12
    't(v;12p)',
    #t15
    't(15;17)(q24.1;q21.2)',
    #t16
    't(16;16)(p13.1;q22)', 't(16;21)(q24.3;q22.1)', 't(16;21)(p11.2;q22.2)',
    #t17
    't(v;17p)', 't(v;17q21.2)',
    #der translocations
    'der(5)t(v;5q)', 'der(7)t(v;7q)',
    'der(12)t(v;12p)', 'der(17)t(v;17p)',
    ]
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
    fish_results=None
    #report = "  44,X,-Y,-Y[17]/46,XY[3]   "
    #report = "47,XY,+21c[6]/48,sl,+11,der(19)t(1;19)(q23;p13.3)[4]"
    #report = "47,XY,+11,t(3;19)(q26.2;p13.3)[4]"
    #report = " 46,xx,t(8;16)(p11.2;p13.3)[20]"
    #report = "45,XX,t(3;21)(q26;q?11.2),del(5)(q23-31q33),-7[14]"
    #report = "46,XX,t(3;17)(q26.2;p13)[18]"
    #report = "46,XY,inv(16)(p13.2q22)[2]/47,sl,del(6)(q13q23),add(5)(q),+22[6]/48,sdl1,+13[4]/46,XY[4]"
    #report = "45,XX,-7[22]/46,idem,+12[3]/47,idem,+12,+20[5]"
    #report = "47,XY,+13,i(13)(q10)x2[2]/47,XY,+13[4]/46,XY[8]"
    #report = "46,XY,del(3)(p13),del(5)(q15q33),del(7)(p13p22),add(12)(p13)[3]/45,sl,dic(20;21)(q1;p1)[5]/47,sl1,+add(21)(p1),+mar[2]"
    #report = "46,XY,t(12;20)(q15;q11.2)[6]/47,sl,+13[2]/94,sdl1x2[2]/93,sdl2,dic(5;6)(q1?1.2;q12~13)[3]/95,sdl2,+15[2]"
    #report = "46,XX,t(12;20),-17,+mar[10]/92,slx2,-t(12;20),-mar[10]"
    #report = "46,XX,t(3;3)(q21.4;q26),inv(3)(q21q26)[20],inv(16)(p13q22.3)"
    #report = "46,XY,t(5;11)(q35;p11)?c,?add(16)(q23~q24)[10]"
    #report = "45,XX,add(1)(p11),-3,add(5)(q31),add(8)(p11),?add(9)(q34),-12,-13,-17,?add(19)(q13),-22,+4mar,inc[cp5]/46,XX[2]"
    #report = "49~51,XY,der(5)t(5;6)(q23;q13),-6,i(9)(q10),+11,del(12)(p12),+19,+22,+2~4mar,inc[21]"
    #report = "47,XY,+21c[6]/48,sl,+11,der(19)t(1;19)(q23;p13.3)[4]"
    #report = "46,XY,t(9;22;10)(q34;q11.2;q11)[12]/45,X,-Y,t(9;22;10)(q34;q11.2;q11)[8]"
    #report = "45,XY,-7,del(13)(q14q31),der(16)t(7;16)(q11.2;q12-13)[4]/46,XY[6]"
    #report = "46,XX,t(9;11;5;22)(p21;q23.3;q35;q12)[10]"
    report = "43~45,XX,-2,?del(3)(q2?),-4,del(4)(p1?),-6,-7,-10,add(12)(p1?),-15,-17,-19,+r,+5~6mar[cp10]"
    result = extract_from_string(report, props, fish=fish_results, verbose = VERBOSE,
                                 only_positive= True)
    print(report)
    print(result)
    #print(props)