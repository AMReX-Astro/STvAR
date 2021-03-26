def parse_name_rank(string):
    if string.find('_') == -1:
        symbname = string
        rank = ''
    else:
        symbname = string[0:string.find('_')] 
        rank = string[string.find('_')+1:]
    return symbname, rank

def unparse_name_rank(symbname, rank = ''):
    if len(rank) != 0:
        string = symbname + '_' + rank
    else:
        string = symbname
    return string

def parse_name_rank_component(string):
    if string.find('_') == -1:
        symbname = string
        rank = ''
        component = ''
    else:
        symbname = string[0:string.find('_')] 
        rank = string[string.find('_')+1:string.rfind('_')]
        component = string[string.rfind('_')+1:]
    return symbname, rank, component

def unparse_name_rank_component(symbname, rank = '', component = ''):
    if len(rank) != 0:
        string = symbname + '_' + rank + '_' + component
    else:
        string = symbname
    return string
    