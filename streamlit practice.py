#!/usr/bin/env python
# coding: utf-8

# In[4]:


import streamlit as st
import pandas as pd
import base64


# In[2]:


st.write("Hello world")


# In[3]:


file = st.file_uploader("upload excel here")


# In[5]:


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


# Examples
df = pd.DataFrame({'x': list(range(10)), 'y': list(range(10))})
st.write(df)

if st.button('Download Dataframe as CSV'):
    tmp_download_link = download_link(df, 'YOUR_DF.csv', 'Click here to download your data!')
    st.markdown(tmp_download_link, unsafe_allow_html=True)

s = st.text_input('Enter text here')
st.write(s)

if st.button('Download input as a text file'):
    tmp_download_link = download_link(s, 'YOUR_INPUT.txt', 'Click here to download your text!')
    st.markdown(tmp_download_link, unsafe_allow_html=True)


# In[ ]:




