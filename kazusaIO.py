from html.parser import HTMLParser
import urllib
import urllib.request

class KazusaHTMLParser(HTMLParser):
    '''
    HTML parser to extract the raw codon usage table from Kazusa HTML result
    '''       
    def __init__(self):
        HTMLParser.__init__(self)
        self.CUTableIncoming = False
        self.res = ""
   
    def handle_starttag(self, tag, attrs):
        if (tag=="pre"): # table is in a <pre>-tag
            self.CUTableIncoming = True
            
    def handle_data(self, data):
        if self.CUTableIncoming: # save content of <pre>
            self.res = data
            self.CUTableIncoming = False
                        
    def getResult(self):
        '''
        return parsed CU table, only gives results != "" after parsing
        '''
        return self.res
    
def getCU(taxid):
    '''
    get codon usage in the form:
    [codon][amino acid][relative frequency of codon] ...
    one line per codon
    '''
    # construct URL for CU of species with given taxid
    req = "http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species="+taxid+"&aa=1"
    res = urllib.request.urlopen(req).read().decode("utf-8") # get HTML
    
    # parse HTML
    p = KazusaHTMLParser()    
    p.feed(str(res))
    
    rstr = p.getResult().replace(")  ", ")\n") # one codon per line
    #remove empty lines
    rlines = rstr.split("\n")
    rlines2 = list()
    for i in range(0, len(rlines)):
        if rlines[i] != "":
            rlines2.append(rlines[i])
    return "\n".join(rlines2)
    