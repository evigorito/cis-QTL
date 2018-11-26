import requests
import lxml.html as lh
import pandas as pd

def get_html(url, out):
    """ Given a url to a html table reads it as pandas data frame and save it as csv file"""
    req=requests.get(url, headers={'User-Agent': 'Mozilla/5.0'})
    doc = lh.fromstring(req.content)
    #Parse data that are stored between <tr>..</tr> of HTML
    tr_elements = doc.xpath('//tr')
    #Create empty list
    col=[]
    i=0
    #For each row, store each first element (header) and an empty list
    for t in tr_elements[0]:
        i+=1
        name=t.text_content()
        #print('%d:"%s"'%(i,name))
        col.append((name,[]))
        ##Since out first row is the header, data is stored on the second row onwards
        #print(len(tr_elements))
    for j in range(1,len(tr_elements)):
        #T is our j'th row
        T=tr_elements[j]
        i=0    
        #Iterate through each element of the row
        for t in T.iterchildren():
            data=t.text_content() 
            #Check if row is empty
            if i>0:
                #Convert any numerical value to integers
                try:
                    data=int(data)
                except:
                    pass
            #Append the data to the empty list of the i'th column
            col[i][1].append(data)
            #Increment i for the next column
            i+=1
    Dict={title:column for (title,column) in col}
    df=pd.DataFrame(Dict)
    df.to_csv(out)

                    
get_html(url=str(snakemake.params['table']), out=str(snakemake.output))
