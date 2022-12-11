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

fetchHosts = {
    "pdb": "http://ftp.wwpdb.org/pub/pdb",
    "pdbe": "ftp://ftp.ebi.ac.uk/pub/databases/pdb",
    "pdbj": "ftp://ftp.pdbj.org/pub/pdb",
}

hostPaths2 = {
    "pdb"  : "https://api.esmatlas.com/fetchPredictedStructure/{code}.pdb",
    "cif"  : "https://api.esmatlas.com/fetchPredictedStructure/{code}.cif",
}

def _fetch3(code, name, state, finish, discrete, multiplex, zoom, type, path, file, quiet, _self=cmd):
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
    nameFmt = '{code}.{type}'
    if type == 'pdb':
        pass
    elif type in ('fofc', '2fofc'):
        nameFmt = '{code}_{type}.ccp4'
    elif type == 'emd':
        nameFmt = '{type}_{code}.ccp4'
    elif type in ('cid', 'sid'):
        bioType = 'pubchem'
        nameFmt = '{type}_{code}.sdf'
    elif type == 'cif':
        pass
    elif type == 'mmtf':
        pass
    elif type == 'cc':
        nameFmt = '{code}.cif'
    elif re.match(r'pdb\d+$', type):
        bioType = 'bio'
    else:
        raise ValueError('type')

    url = hostPaths2[bioType]
    url_list = []
    for url in url if cmd.is_sequence(url) else [url]:
        url_list += [url] if '://' in url else [fetch_host + url for fetch_host in fetch_host_list]

#    if bioType not in ['cc','bcif']:
#        code = code.lower()

    fobj = None
    contents = None

    if not file or file in (1, '1', 'auto'):
        file = os.path.join(path, nameFmt.format(code=code, type=type))

    if not is_string(file):
        fobj = file
        file = None
    elif os.path.exists(file):
        # skip downloading
        url_list = []

    for url in url_list:
        url = url.format(mid=code[-3:-1], code=code, type=type)

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

def _multifetch3(code,name,state,finish,discrete,multiplex,zoom,type,path,file,quiet,_self):
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

        r = _fetch3(obj_code, obj_name, state, finish,
                   discrete, multiplex, zoom, type, path, file, quiet, _self)

        if chain and isinstance(r, str):
            if _self.count_atoms(r'?%s & c. \%s' % (r, chain)) == 0:
                _self.delete(r)
                raise pymol.CmdException('no such chain: ' + chain)
            _self.remove(r'?%s & ! c. \%s' % (r, chain))

    return r

import urllib.request
from bs4 import BeautifulSoup     
def esmatlas_fetch(code, name='', state=0, finish=1, discrete=-1,
           multiplex=-2, zoom=-1, type='pdb', async_=0, path='',
           file=None, quiet=1, _self=cmd, **kwargs):

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
#  fp = urllib.request.urlopen("https://www.uniprot.org/uniprot/"+code) # gone outdated with uniprot updates
#https://rest.uniprot.org/uniprotkb/P19484.txt somtimes toolarge to fetch
#https://rest.uniprot.org/uniprotkb/P19484.fasta reasonable and informative
    #fp = urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/"+code+".fasta")
    #mybytes = fp.read()
    #mystr = mybytes.decode("utf8")
    #fp.close()
    #html=mystr
    #soup = BeautifulSoup(html,features="html.parser") 
    #print(soup.getText())
   
    if async_:
        _self.async_(_multifetch3, *args, **kwargs)
    else:
        try:
            _self.block_flush(_self)
            r = _multifetch3(*args)
        finally:
            _self.unblock_flush(_self)
    return r

pymol.cmd.extend("esm", esmatlas_fetch)
