#!/usr/bin/env python
# coding: utf-8

# In[2]:


import re
import pandas as pd
import streamlit as st
import base64


# # identifying karyotype reports

# In[ ]:


st.write('# Cytogenetics calculator')


# In[2]:


#karyotypes = pd.read_excel('Cytogenetics_TS_Apr2021.xlsx', )
file = st.file_uploader("Upload excel here")
if file:
    karyotypes = pd.read_excel(file)
    karyotypes = karyotypes.drop(columns='ID')
    #karyotypes.head()


    # In[4]:


    properties = karyotypes.columns[4:-1].to_list()


    # In[5]:


    def properties_dict(properties=properties):
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
                f = v[:m.start()] + '(' + m.group()[:-1] +                 ')(' + m.group()[-1]+ v[m.end():]
                v=f
                
            #creating escape characters for strings
            v = re.escape(v)
            
            if v == 't\\(v;11\\)':
                v = 't\\(\d+;11\\)'
            
            d[k] = v
        return {v:k for k,v in d.items()}


    # In[6]:


    prop_dict = properties_dict()
    #prop_dict


    # In[7]:


    def remove_artefact(row):
        if re.search('_x000D_$', row['Cytogenetics']):
            return row['Cytogenetics'].strip()[:-7]
        return row['Cytogenetics'].strip()


    # In[8]:


    karyotypes['Cytogenetics'] = karyotypes.apply(remove_artefact, axis=1)


    # In[9]:


    karyotypes['Error'] = None


    # In[10]:


    def gram_error(string):
        punct = r'()[]/'
        punct_dict = {s:0 for s in punct}
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
        pass


    # In[11]:


    def parse_karyotype(row):
        error = gram_error(row['Cytogenetics'])
        if error:
            row['Error'] = error
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
        row['Number of cytogenetic abnormalities'] = len(abnorms)
        row['Monosomy'] = mono
        row['Structural'] = struc
        row['abnormal(17p)'] = seventeen_p
        for c in col_true:
            row[prop_dict[c]] = True
        return row


    # In[12]:


    results = karyotypes.apply(parse_karyotype, axis=1)
    #results.head()


    # In[13]:


    #results.any()


    # In[15]:


    #results.to_excel('Cytogenetics_output_V2.xlsx')


    # In[ ]:


    def download_link(object_to_download, download_filename, download_link_text):
        """
        Generates a link to download the given object_to_download.

        object_to_download (str, pd.DataFrame):  The object to be downloaded.
        download_filename (str): filename and extension of file. e.g. mydata.csv, some_txt_output.txt
        download_link_text (str): Text to display for download link.

        Examples:
        download_link(YOUR_DF, 'YOUR_DF.csv', 'Click here to download data!')
        download_link(YOUR_STRING, 'YOUR_STRING.txt', 'Click here to download your text!')

        """
        if isinstance(object_to_download,pd.DataFrame):
            object_to_download = object_to_download.to_csv(index=False)

        # some strings <-> bytes conversions necessary here
        b64 = base64.b64encode(object_to_download.encode()).decode()

        return f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>'


    # In[ ]:


    st.write(results)


    # In[ ]:


    if st.button('Download Dataframe as CSV'):
        tmp_download_link = download_link(results, 'YOUR_DF.csv', 'Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)

