import json

class JsonFile:
    '''Compare a json object againt the specified schema.
    '''
    old_schema = {
            'submitter': {
                'team':''
                ,'name':''
            }
            ,'vcf':''
            ,'geneList':[{'gestalt_score':''}]
            ,'features':''
            ,'ranks':''
            ,'genomicData':[
                {   'Mutations':{'HGVS-code':'','Mutation 1':''}
                    ,'Test Information':{'Mutation Type':''}
                }
            ]
    }

    schema = {'algo_deploy_version': '',
            'case_id': '',
            'detected_syndromes': [{'combined_score': '',
                'feature_score': '',
                'gestalt_score': '',
                'has_mask': '',
                'omim_id': '',
                'syndrome_name': ''}],
            'documents': [],
            'features': [],
            'genomic_entries': [],
            'selected_syndromes': [{'has_mask': '', 'omim_id': '', 'syndrome_name': ''}],
            'submitter': {'user_email': '', 'user_name': '', 'user_team': ''}}


    def __init__(self, path=None, inputdict=None):
        self.process = {
            'genomic_entries': [self.load_entries]
            }
        self.raw = None
        self.errors = None
        if path is not None:
            self.read_file(path)
        elif inputdict is not None:
            self.read_dict(inputdict)

    def load_entries(self, filename, path):
        return json.load(open(os.path.join(path,filename),'r'))

    def process(self):
        self.check()
        self.raw

    def read_file(self, inputpath):
        read_json = json.load(open(inputpath,'r'))
        self.raw = read_json

    def read_dict(self, inputdict):
        self.raw = inputdict

    def generate(self):
        return self.generate_schema(self.raw)

    def generate_schema(self, data):
        '''Recursively get schema from available raw data'''
        if isinstance(data, dict):
            return {k:self.generate_schema(v) for k,v in data.items()}
        elif isinstance(data, list):
            res = [self.generate_schema(v) for v in data]
            out = []
            for r in res:
                if isinstance(r, dict):
                    if set(r.keys()) not in [set(o) for o in out]:
                        out.append(r)
                elif isinstance(r, str):
                    if r not in out:
                        out.append(r)
                else:
                    out.append(r)
            return out
        else:
            return ''

    def filter_schema(self, schema):
        if isinstance(schema,dict):
            out = {}
            for k,v in schema.items():
                if isinstance(v, tuple):
                    if not v[0]:
                        out[k] = v
                else:
                    o = self.filter_schema(v)
                    if not(o == {} or o == []):
                        out[k] = o
            return out
        elif isinstance(schema,list):
            out = []
            for v in schema:
                if isinstance(v,tuple):
                    if not v[0]:
                        out.append(v)
                else:
                    o = self.filter_schema(v)
                    if not(o =={} or o ==[]):
                        out.append(o)
            return out
        else:
            print(schema)
            raise RuntimeError()


    def check(self, false_only=True):
        schema =  self.check_schema(self.schema, self.raw)
        self.error = schema
        if false_only:
            return self.filter_schema(schema)
        else:
            return schema

    def check_schema(self, schema, rawjs):
        '''Recursively check the specified schema against the provided schema
        '''
        if isinstance(schema, dict):
            if not isinstance(rawjs,dict):
                return (False, "Not a dict")
            else:
                out = {}
                for k,v in schema.items():
                    if k not in rawjs:
                        out[k] = (False, "No key")
                    else:
                        out[k] = self.check_schema(v,rawjs[k])
                return out
        elif isinstance(schema, list):
            if not isinstance(rawjs,list):
                return (False, "Not a list")
            else:
                if len(schema) == 0:
                    return (True, "")
                out = []
                for v in rawjs:
                    res = []
                    for s in schema:
                        res.append(self.check_schema(s,v))
                    if len(schema) == 1:
                        res = res[0]
                    out.append(res)
                return out
        elif hasattr(schema, '__call__'):
            return schema(rawjs)
        else:
            if rawjs is None:
                return (False, "No value")
            elif schema == '':
                return (True, "")
            elif schema == rawjs:
                return (True, "matches expected string")
            else:
                return (False,"No value")
