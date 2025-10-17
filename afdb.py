import re
import os
import sys
import copy
import traceback
import pymol
cmd = sys.modules["pymol.cmd"]
from . import selector
from . import colorprinting
from .cmd import _cmd, \
      DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, \
      is_list, space_sc, safe_list_eval, is_string, loadable
from .constants import _loadable
from pymol.creating import unquote


import re

def parse_id(s):
    if "AF_AF" == s[0:5] and "F1" == s[-2:]:
        # RCSB's AFDB copy like "AF_AFA0A009IHW8F1"
        print("seems RCSB")
        return (s[5:-2])
    if "AF-" == s[0:3]:
        # AFDB ID, especially model file names like "AF-Q5VSL9-F1-model_v6.cif"
        print("seems AFDB ID")
        ss = s.split("-")
        if ss[0] == "AF" and ss[2] == "F1" and ss[3].split("_")[0] == "model":
            return ss[1]
    else:
        if "-" in s:
            ss = s.split("-")
            for sss in ss:
                if "." in sss:
                    sss = sss.split(".")[0]
                if (len(sss) == 5 or len(sss) == 10):
                    cs = re.sub(r'[^0-9]', '', sss)
                    if all(c.isupper() for c in cs):
                        return(sss)
                    
        else:
            clean = re.sub(r'[^A-Z0-9]', '', s)
            return clean

def extract_uniprotID(s):
    if "://" in s:
        # seems like URL
        ss=s.split("/")
        for s in ss:
            # UniProt ID has Upper Cases, but other URL elements would not
            if any(c.isupper() for c in s):
                return parse_id(s)
    else:
        return parse_id(s)

def extract_middle_token(s):
    # Split on '-' or '/'
    parts = re.split(r'[-/]', s)
    # Return the middle part if there are 3 or more segments
    if len(parts) >= 3:
        return parts[len(parts)//2]
    # Otherwise, just return the last meaningful part
    return parts[-1]

fetchHosts = {
    "pdb": "http://ftp.wwpdb.org/pub/pdb",
    "pdbe": "ftp://ftp.ebi.ac.uk/pub/databases/pdb",
    "pdbj": "ftp://ftp.pdbj.org/pub/pdb",
}

hostPaths = {
    "pdb"  : "https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v{version}.pdb",
    "cif"  : "https://alphafold.ebi.ac.uk/files/AF-{code}-F1-model_v{version}.cif",
}

versions = [6,5,4,3,2,1]

def _fetch2(code, name, state, finish, discrete, multiplex, zoom, type, path, file, quiet, _self=cmd):
    '''
        code = str: single pdb identifier
        name = str: object name
        state = int: object state
        finish =
        discrete = bool: make discrete multi-state object
        multiplex = bool: split states into objects (like split_states)
        zoom = int: zoom to new loaded object
        type = str: fofc, 2fofc, pdb, pdb1, ... 
        path = str: fetch_path
        file = str or file: file name or open file handle
        '''
    r = DEFAULT_ERROR
    fetch_host_list = [x if '://' in x else fetchHosts[x]
                       for x in _self.get("fetch_host").split()]

    # file types can be: fofc, 2fofc, pdb, pdb1, pdb2, pdb3, etc...
    # bioType is the string representation of the type
    # nameFmt is the file name pattern after download
    bioType = type
    nameFmt = 'AF-{code}-F1-model_v{version}.{type}'
    url = hostPaths[bioType]
    url_list = []
    for url in url if cmd.is_sequence(url) else [url]:
        for version in versions:
            format_url = url.format(mid=code[-3:-1], code=code, type=type, version=str(version))
            url_list += [format_url] if '://' in url else [fetch_host + url for fetch_host in fetch_host_list]

#    if bioType not in ['cc','bcif']:
#        code = code.lower()

    fobj = None
    contents = None
    
    for url in url_list:
        print(url)
        filename = url.split("/")[-1]        
        if name is not None:
            name = filename.split(".")[0]
        file = os.path.join(path, filename)

        if not is_string(file):
            fobj = file
            file = None
        elif os.path.exists(file):
            # skip downloading
            pass            
        try:
            contents = _self.file_read(url)

            # assume HTML content means error on server side without error HTTP code
            if b'<html' in contents[:500].lower():
                raise pymol.CmdException

        except pymol.CmdException:
            if not quiet:
                colorprinting.warning(" Warning: failed to fetch from %s" % (url,))
            continue

        if file:
            try:
                fobj = open(file, 'wb')
            except IOError:
                colorprinting.warning(' Warning: Cannot write to "%s"' % file)

        if fobj:
            fobj.write(contents)
            fobj.flush()
            if file:
                fobj.close()

        if not file:
            return DEFAULT_SUCCESS

        break

    if os.path.exists(file):
        r = _self.load(file, name, state, '',
                       finish, discrete, quiet, multiplex, zoom)
    elif contents and bioType in ('pdb', 'bio'):
        r = _self.read_pdbstr(contents, name, state,
                              finish, discrete, quiet, zoom, multiplex)
    elif contents and bioType in ('cif', 'cc'):
        r = _self.load_raw(contents, 'cif', name, state,
                           finish, discrete, quiet, multiplex, zoom)
    elif contents and bioType in ('mmtf',):
        r = _self.load_raw(contents, 'mmtf', name, state,
                           finish, discrete, quiet, multiplex, zoom)

    if not _self.is_error(r):
        return name

    colorprinting.error(" Error-fetch: unable to load '%s'." % code)
    return DEFAULT_ERROR

def _multifetch2(code,name,state,finish,discrete,multiplex,zoom,type,path,file,quiet,_self):
    import string
    r = DEFAULT_SUCCESS
    code_list = code.split()
    name = name.strip()
    if (name!='') and (len(code_list)>1) and (discrete<0):
        discrete = 1 # by default, select discrete  when loading
        # multiple PDB entries into a single object

    all_type = type
    for obj_code in code_list:
        obj_name = name
        type = all_type

        if not type:
            if 1 < len(obj_code) < 4:
                type = 'cc'
            else:
                type = _self.get('fetch_type_default')

        # allow fetching codes like EMD-3489 or emd_3489
        #if obj_code[3:4] in ('_', '-') and \
        #   obj_code[:3].upper() in ('CID', 'SID', 'EMD'):
        #    if not obj_name:
        #        obj_name = obj_code
        #        type = obj_code[:3].lower()
        #        obj_code = obj_code[4:]

        if not obj_name:
            obj_name = obj_code
            if type.endswith('fofc'):
                obj_name += '_' + type
            elif type == 'emd':
                obj_name = 'emd_' + obj_code

        chain = None
        #if len(obj_code) in (5,6,7) and type in ('pdb', 'cif', 'mmtf'):
        #    obj_code = (
        #        obj_code
        #        .replace('.', '')
        #        .replace('_', '')
        #        .replace('-', '')
        #        .replace(':', '')
        #    )
        #    obj_code, chain = obj_code[:4], obj_code[4:]

        obj_name = _self.get_legal_name(obj_name)

        r = _fetch2(obj_code, obj_name, state, finish,
                   discrete, multiplex, zoom, type, path, file, quiet, _self)

        if chain and isinstance(r, str):
            if _self.count_atoms(r'?%s & c. \%s' % (r, chain)) == 0:
                _self.delete(r)
                raise pymol.CmdException('no such chain: ' + chain)
            _self.remove(r'?%s & ! c. \%s' % (r, chain))

    return r

import urllib.request
from bs4 import BeautifulSoup     
def afdb_fetch(code, name='', state=0, finish=1, discrete=-1,
           multiplex=-2, zoom=-1, type='cif', async_=0, path='',
           file=None, quiet=1, _self=cmd, **kwargs):
    code = extract_uniprotID(code)
    state, finish, discrete = int(state), int(finish), int(discrete)
    multiplex, zoom = int(multiplex), int(zoom)
    async_, quiet = int(kwargs.pop('async', async_)), int(quiet)

    if kwargs:
        raise pymol.CmdException('unknown argument: ' + ', '.join(kwargs))

    r = DEFAULT_SUCCESS
    if not path:
        # blank paths need to be reset to '.'
        path = _self.get('fetch_path') or '.'
    if async_ < 0: # by default, run asynch when interactive, sync when not
        async_ = not quiet
    args = (code, name, state, finish, discrete, multiplex, zoom, type, path, file, quiet, _self)
    kwargs = { '_self' : _self }


    # Print Unipot paget title
    uni = "https://rest.uniprot.org/uniprotkb/"+code+".txt"
    with urllib.request.urlopen(uni) as fp:
        data = fp.read().decode("utf-8")
        # Print content
        print(data)

    # Save to file named "<code>.txt"
    if len(data) == 0:
        print(f"This entry {code} is not available in UniProt")
        uni = "https://rest.uniprot.org/unisave/" + code
        with urllib.request.urlopen(uni) as fp:
            data = fp.read().decode("utf-8")
            print(data)

    with open(f"{code}.txt", "w", encoding="utf-8") as f:
        f.write(data)
        print(f"Saved as {code}.txt")
   
    if async_:
        _self.async_(_multifetch2, *args, **kwargs)
    else:
        try:
            _self.block_flush(_self)
            r = _multifetch2(*args)
        finally:
            _self.unblock_flush(_self)
    return r

pymol.cmd.extend("af", afdb_fetch)
