#!/usr/bin/env python
# coding: utf-8

# # Sample Jupyter Notebook
# 
# This notebook checks that pMuTT was installed correctly on your computer. 
# 
# The following checks that you have the most updated version. Click the cell below and click the **run** button. The expected result is:
# 1.2.11

# In[ ]:


import pmutt

print(pmutt.__version__)

