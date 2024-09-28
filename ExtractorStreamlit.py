#importing modules

import pandas as pd
import streamlit as st
from Extractor import *



# def download_link(object_to_download, download_filename, download_link_text):
#         """
#         Generates a link to download the given object_to_download.

#         object_to_download (str, pd.DataFrame):  The object to be downloaded.
#         download_filename (str): filename and extension of file. e.g. mydata.csv, some_txt_output.txt
#         download_link_text (str): Text to display for download link.

#         Examples:
#         download_link(YOUR_DF, 'YOUR_DF.csv', 'Click here to download data!')
#         download_link(YOUR_STRING, 'YOUR_STRING.txt', 'Click here to download your text!')

#         """
#         if isinstance(object_to_download,pd.DataFrame):
#             object_to_download = object_to_download.to_csv(index=False)

#         # some strings <-> bytes conversions necessary here
#         b64 = base64.b64encode(object_to_download.encode()).decode()

#         return f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>'


st.write('# Cytogenetics calculator')

base_mode = st.radio('For the latest version, select \'ELN2022\'', ['ELN2022','BJH2021'])

prop_dict = available_configs()[base_mode]

st.write('## text input')

cytogenetics = st.text_input('type/paste cytogenetic report here')
if cytogenetics:
    st.write(cytogenetics)
    result = extract_from_string(cytogenetics, prop_dict, only_positive=True)
    st.write('#### Results')
    if result['error']:
        st.write('ERROR PRESENT:')
        st.write('-'*20)
        for m in result['error_message']:
            st.write(m)
    else:
        st.write('-'*20)
        if len(result['Warnings'])>0:
            st.write('### WARNING')
            for w in result['Warnings']:
                st.write(w)
            st.write('-'*20)
            
        st.write('Number of cytogenetic abnormalities: ', 
                 result['result']['Number of cytogenetic abnormalities'])
        st.write('Polysomy count:', 
                 result['result']['Polysomy'])
        st.write('Monosomy count:', 
                 result['result']['Monosomy'])
        st.write('Structural abnormality count: ', 
                 result['result']['Structural'])
        
        non_count_results = [r for r in result['result'] if r not in 
                             ['Number of cytogenetic abnormalities',
                              'Polysomy','Monosomy', 'NonSexChromosomeMonosomies',
                              'Structural', 'Warnings',
                              'Cytogenetics', 'Error', 'Error description']]
        if len(non_count_results)>0:    
            st.write(f'Other abnormalities detected from pre-set list: {non_count_results}')
        st.write('-'*20)
        
        
        

# st.write('## file input')
# file = st.file_uploader("Upload CSV file here")
# if file:
#     karyotypes = load_file(file)
#     karyotypes['Error'] = bool(False)
#     karyotypes['Error description'] = None
#     results = karyotypes.apply(parse_karyotype_clone, axis=1, args= (prop_dict,))
#     results.loc[results['Error'] == False] = results.loc[results['Error']==False].fillna(False)
#     results['Error'] = results['Error'].astype(bool)
    
#     #sorting column order
#     first_order_cols = ['Cytogenetics', 'Error', 'Error description', 'Warnings',
#              'Number of cytogenetic abnormalities', 'Polysomy','Monosomy',
#                         'NonSexChromosomeMonosomies', 'Structural']
#     rest_cols = [c for c in results.columns if c not in first_order_cols]
#     results = results[first_order_cols + rest_cols]
    #st.dataframe(results)
    
    #st.download_button(
    #     label="Click to download results as CSV",
    #     data=results
    #     file_name='cytogenetics_results.csv',
    #     mime='text/csv',
    # )
    #if st.button('Download Dataframe as CSV'):
    #    tmp_download_link = download_link(results, 'YOUR_DF.csv', 'Click here to download your data!')
    #    st.markdown(tmp_download_link, unsafe_allow_html=True)
st.info('Uploaded data is never saved')
