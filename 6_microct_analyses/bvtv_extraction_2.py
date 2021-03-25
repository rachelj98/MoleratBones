# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:52:53 2018
Note that this script must be in same directory as files to work. -RJ
@author: annal
"""


import PyPDF2
import xlsxwriter
import os

def clean(line):
    l = line.strip()
    for c in l:
        if not c.isnumeric() and not c=='.' and not c=='-':
            lastInd = l.find(c)
            break
    
    return float(l[:lastInd])

def getInfo(page):
    lines = page.split('\n')
    lines = lines[12:]
    direct_data = {'TV':clean(lines[2]),
                   'BV':clean(lines[3]),
                   'BV/TV':clean(lines[4]),
                   'Conn. D.':clean(lines[5]),
                   'SMI':clean(lines[6]),
                   'TB.N':clean(lines[16]),
                   'TB.Th':clean(lines[17]),
                   'TB.Sp':clean(lines[18])
                   }
    TRI_data = {'TV':clean(lines[7]),
                   'BV':clean(lines[8]),
                   'BV/TV':clean(lines[9]),
                   'BS':clean(lines[10]),
                   'BS/BV':clean(lines[11]),
                   'TB.N':clean(lines[19]),
                   'TB.Th':clean(lines[20]),
                   'TB.Sp':clean(lines[21])
                   }
    anisotropy = {'H1':clean(lines[12]),
                   'H2':clean(lines[13]),
                   'H3':clean(lines[14]),
                   'DA':clean(lines[15])}
    return direct_data, TRI_data, anisotropy

def extract(pdf_name):
    pdfFileOBJ = open(pdf_name,'rb')
    pdfReader = PyPDF2.PdfFileReader(pdfFileOBJ)
    pageObj = pdfReader.getPage(0)
    page = pageObj.extractText()
    data = {'direct':getInfo(page)[0],'TRI':getInfo(page)[1],'an':getInfo(page)[2]}
    pdfFileOBJ.close()
    return data
    
pdfs = []
worksheets = []
for file in os.listdir('/directory_to_folder_containing_PDFs_provided_by_vivaCTscanner'):
    if file[-3:]=="PDF":
        pdfs += [file]
workbook = xlsxwriter.Workbook('bvtv_lv6_v2_1dec2019.xlsx')
worksheet = workbook.add_worksheet()
pos = [0,0]
headers = ['Bone']
data = extract(pdfs[0])
bone = pdfs[0].split('_')[1]
for sub in data:
        d = data.get(sub)
        for k in d:
            headers += [k]
for header in headers:
    worksheet.write(pos[0],pos[1],header)
    pos[1] += 1
        
pos = [1,0]
for pdf in pdfs:
    data = extract(pdf)
    numbers = [pdf.split('_')[1] + pdfs[0].split('_')[2]]
    for sub in data:
        d = data.get(sub)
        for k in d:
            numbers += [d.get(k)]
    for value in numbers:
        worksheet.write(pos[0],pos[1],value)
        pos[1] += 1
    pos[0] += 1
    pos[1] = 0
        
        
workbook.close()
    