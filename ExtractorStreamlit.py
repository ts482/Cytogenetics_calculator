#importing modules

import re
import pandas as pd
import streamlit as st
import base64
from Extractor import *

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


st.write('# Cytogenetics calculator')
file = st.file_uploader("Upload CSV file here")
if file:
    karyotypes = load_file(file)
    prop_dict = properties_dict(karyotypes=karyotypes)
    karyotypes['Cytogenetics'] = karyotypes.apply(remove_artefact, axis=1)
    karyotypes['Error'] = bool(False)
    karyotypes['Error description'] = None
    results = karyotypes.apply(parse_karyotype, axis=1, args= (prop_dict,))
    results.loc[results['Error'] == False] = results.loc[results['Error']==False].fillna(False)
    results['Error'] = results['Error'].astype(bool)
    st.dataframe(results)
    
    if st.button('Download Dataframe as CSV'):
        tmp_download_link = download_link(results, 'YOUR_DF.csv', 'Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)

st.info('Uploaded data is never saved')