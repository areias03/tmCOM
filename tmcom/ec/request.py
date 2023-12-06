def modelseed_query(seed_id: str):
    '''
    Queries the ModelSEED database for other id's for BIGG and KEGG.

    :param seed_id: A ModelSEED entry ID.
    :return: Entry ID for BIGG, list of entry ID's for KEGG.
    :rtype: str,str
    '''
    SOLR_URL='https://modelseed.org'
    
    if seed_id == None:
        bigg_id = 'None'
        kegg_id = 'None'
    else:
        try:
            connection = urlopen(SOLR_URL+f'/solr/reactions/select?wt=json&q=id:{seed_id}&fl=name,id,formula,charge,aliases')
            response = json.load(connection)
            for document in response['response']['docs']:
                bigg_id = list(filter(lambda a: 'BiGG:' in a, document.get('aliases')))
                kegg_id = list(filter(lambda a: 'KEGG:' in a, document.get('aliases')))
                if len(bigg_id)== 0 and len(kegg_id)== 0:
                    bigg_id = 'None'
                    kegg_id = 'None'
                elif len(bigg_id)== 0 and len(kegg_id)!= 0:
                    bigg_id = 'None'
                    kegg_id = list(kegg_id)[0]
                    kegg_id = kegg_id.replace('KEGG: ','')
                elif len(bigg_id)!= 0 and len(kegg_id)== 0:
                    kegg_id = 'None'
                    bigg_id = list(bigg_id)[0]
                    bigg_id = bigg_id.replace('BiGG: ','')
                else:
                    kegg_id = list(kegg_id)[0]
                    kegg_id = kegg_id.replace('KEGG: ','')
                    bigg_id = list(bigg_id)[0]
                    bigg_id = bigg_id.replace('BiGG: ','')
        except:
            bigg_id = 'None'
            kegg_id = 'None'

    return bigg_id,kegg_id

def bigg_query(bigg: str):
    '''
    Queries the BIGG database for EC numbers.
    
    :param bigg: A BIGG entry ID.
    :return: List of EC numbers.
    :rtype: list
    '''
    
    bigg_ls = []
    url =f'http://bigg.ucsd.edu/api/v2/universal/reactions/{bigg}'

    with requests.request("GET", url) as resp:
        try:
            resp.raise_for_status()  # raises exception when not a 2xx response
            if resp.status_code != 204:
                data = dict(resp.json())
                ec_l = data['database_links']
                if ec_l == None:
                    bigg_ls.append(None)
                else:
                    ec = [i['id'] for i in ec_l['EC Number']]
                    if ec == None:
                        bigg_ls.append(None)
                    bigg_ls.append(ec)
            else: 
                bigg_ls.append(None)
        except:
            bigg_ls.append(None)

    return bigg_ls