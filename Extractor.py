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
        
        #recognising monosomies
        if v.startswith('Monosomy'):
            v = '-' + v[8:]
            
        matches = re.finditer('(\d+)(p|q)',v)
        for m in matches:
            f = v[:m.start()] + '(' + m.group()[:-1] + \
                ')(' + m.group()[-1]+ v[m.end():]
            v=f
            
        #creating escape characters for strings
        v = re.escape(v)
        
        if v == 't\\(v;11\\)':
            v = 't\\(\d+;11\\)'
        
        d[k] = v
    return {v:k for k,v in d.items()}

def remove_artefact(row):
    '''
    Removes whitespaces and import artefacts from reports 
    Function to be used in conjunction with pd.DataFrame.apply, 
    and the column containing the cytogenetic report
    '''
    if re.search('_x000D_$', row['Cytogenetics']):
        return row['Cytogenetics'].strip()[:-7]
    return row['Cytogenetics'].strip()

def gram_error(string):
    '''
    Recognise if there is a grammatical error within a karyotype report
    '''
    missing = list()
    for k in punct_dict:
        punct_dict[k] = len(re.findall(re.escape(k), string))
    if punct_dict['/'] != punct_dict['[']-1:
        return 'incorrrect ratio of \'\\\' to \'[]\' '
    if punct_dict['['] != punct_dict[']']:
        missing.append(min(['[', ']'], key=punct_dict.get))
    if punct_dict['('] != punct_dict[')']:
        missing.append(min(['(',')'], key=punct_dict.get))
    if len(missing)!=0:
        return f'missing grammar: {", ".join(missing)}'
    substring = re.split('/', string)
    for i, s in enumerate(substring):
        chrom = s[:s.index(',')] #chromsome string
        if not re.search('idem', s):
            expected = 46
        expected -= len(re.findall('\-', s))
        expected += len(re.findall('\+', s))
        if re.search('[~]', s):
            low_num, high_num = int(chrom[:2]), int(chrom[-2:]) 
            if low_num <= expected <= high_num:
                pass
            elif expected < low_num:
                return f'chromsome number lower than expected in {i+1} subsection'
            else:
                return f'chromsome number higher than expected in {i+1} subsection'
        else:
            num = int(chrom[:2])
            if expected == num:
                pass
            elif expected > num:
                return f'chromsome number higher than expected in {i+1} subsection'
            else:
                return f'chromsome number lower than expected in {i+1} subsection'
    pass


def parse_karyotype(row, prop_dict):
    '''
    Calculates which cytogenetic abnormalities are present in each report
    and fills in the relevant columns
    '''
    error = gram_error(row['Cytogenetics'])
    if error:
        row['Error'] = True
        row['Error description'] = error
        if not error.startswith('chromsome'):
            return row
    abnorms = set(re.split('/|,|\[', row['Cytogenetics']))
    removed = set()
    col_true = set()
    mono = 0
    struc = 0
    der = 0
    seventeen_p = False
    for a in abnorms:
        if re.fullmatch('\d\d([~-]\d\d)?|X[XY]?|(cp)?\d\d?\]|idem', a):
            removed.add(a)
#         if re.search('mar', a): #no longer the case to remove mar
#             removed.add(a)
        
        if a not in removed:
            
            if re.fullmatch('-\d+|-[XY]', a):
                mono += 1
            if re.search('[inv|t].*\).*\)', a):
                struc += 1
            
            if re.search('der', a):
                der += 1

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
            
            for p in prop_dict:
                if re.search(p, a):
                    col_true.add(p)
    abnorms = abnorms.difference(removed)
    row['Number of cytogenetic abnormalities'] = len(abnorms) + der
    row['Monosomy'] = mono
    row['Structural'] = struc
    row['abnormal(17p)'] = seventeen_p
    for c in col_true:
        row[prop_dict[c]] = True
    return row



if __name__ == '__main__':
    karyotypes = load_file()
    prop_dict = properties_dict(karyotypes=karyotypes)
    karyotypes['Cytogenetics'] = karyotypes.apply(remove_artefact, axis=1)
    karyotypes['Error'] = False
    karyotypes['Error description'] = None
    results = karyotypes.apply(parse_karyotype, axis=1, args=(prop_dict,))
    results.loc[results['Error']==False] = results.loc[results['Error']==False].fillna(False)
    results['Error'] = results['Error'].astype(bool)
    results.to_excel('Cytogenetics_output_V4.xlsx')