'''
Retrieve linear motifs from ELM database.
Created on 15 Mar 2011

@author: Fergal
'''

#Non-built in packages
import SOAPpy
from SOAPpy import WSDL


# web service client configuration   
endpoint = 'http://api.bioinfo.no/services/ELMdb'
namespace = 'http://elm.eu.org/ELMdb'
#wsdl = "http://elm.eu.org/webservice/ELMdb.wsdl"
elm_db = SOAPpy.SOAPProxy(endpoint)
elm_db.namespace = namespace
elm_db.noroot = 1
#server = WSDL.Proxy(wsdl)

def get_all_elms():
    return elm_db.getAllELMs()

def get_elm_regexes():
    return [elm['Regex'] for elm in get_all_elms()]
	
def get_elm_instance(accession):
    return elm_db.getELMInstance(accession)


if __name__ == '__main__':
    print get_elm_regexes()
