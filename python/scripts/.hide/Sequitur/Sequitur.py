'''
DESCRIPTION
    An object which holds all the characters that two entries into the Suffix Tree
    have in common
INPUT
    stalk | a string
METHODS
    common_substring(other) -> (Stalk,Stalk,Stalk)
'''
class Stalk:
    def __init__(self,stalk):
        self.stalk = stalk
        if len(self.stalk) == 0: self.stalk = '$'

    def __repr__(self): return self.stalk

    def __eq__(self,other): return self[0] == other[0]

    def __hash__(self): return hash(self[0])

    def __getitem__(self,index): return self.stalk[index]

    def __len__(self):
        if self.stalk == '$' or self.stalk == '^': return 0
        return len(self.stalk)
    
    def __str__(self):
        if self.stalk =='$': return ''
        else: return self.stalk

    '''
    DESCRIPTION
        A method which finds all the characters two substrings have in common
    INPUT
        other | a Stalk
    OUTPUT
        A 3-tuple of Stalk. The first stalk is what the two stalks have in common,
        the second stalk is what the stalk making the call has proceeding the common
        string, the third stalk is what the other stalk has proceeding the common string
    '''
    def common_substring(self,other):
        i = 0
        substr = ''
        if type(other) == str: other = Stalk(other)
        while i < min(len(self),len(other)) and self[i] == other[i]:
            substr += self[i]
            i += 1
        return Stalk(substr), Stalk(self[i:]), Stalk(other[i:])

'''
DESCRIPTION
INPUT
OUTPUT
'''
class Leaf:
    def __init__(self,left,right=''):
        if len(left) == 0:
            self.left = left
            self.right = 1
        else:
            self.left = left
            self.right = Leaf(right)
    
    def __repr__(self): return str(self.right)

    def __eq__(self,other): return self[0] == other[0]

    def __hash__(self): return hash(self.left)

    def __getitem__(self,index): return self.left[index]

    def __len__(self):
        if self.left == '$': return 0
        return len(self.left)
    
    def __is_shallow__(self): return True

    def reads(self): return set()

'''
DESCRIPTION
INPUT
OUTPUT
    '''
class Branch:
    def __init__(self):
        self.b = {}
        self.s = {}

    def __repr__(self): return repr(self.b)

    def __str__(self):
        s = '' 
        for i in range(len(list(self.b.values()))-1):s+=str(list(self.b.values())[i])+'\n'
        return s+str(list(self.b.values())[-1])
    
    def __getitem__(self,index):
        if type(index) == str: return self.b[Stalk(index)]
        return self.b[index]

    '''
    DESCRIPTION
    INPUT
    OUTPUT
    '''
    def __is_shallow__(self):
        for a in self.b.values():
            if type(a) == Branch: return False
        return True

    '''
    DESCRIPTION
    INPUT
    OUTPUT
    '''
    def __traverse__(self,context):
        b = self[context[0]]
        s = self.s[context[0]]
        context = context[len(s[0]):]
        while len(context) > 0 and len(b) > 1:
            s = b.s[context[0]]
            b = b[context[0]]
            context = context[len(s[0]):]
        return b
    
    def __setitem__(self,index,value):
        if type(index) == str: self.b[Stalk(index)] = value
        else: self.b[index] = value

    def __contains__(self,other): 
        if type(other) == str: return Stalk(other) in self.b
        return other in self.b

    def __len__(self): return len(self.b)

    def pop(self,index): return self.b.pop(index)

    '''
    DESCRIPTION
        adds a suffix to the trie
    INPUT
        stalk | a Stalk() which is a common substring of every read up to this point and beyond
        reads | a set of reads which have with the same common substring up to this point
    '''
    def add(self,stalk,reads):
        if stalk in self:
            if not len(stalk):
                self[stalk].right+=1
                self.s[stalk][1].update(reads)
                return
            if type(self[stalk]) == Leaf:
                branch = Branch()
                l1 = self.pop(stalk)
                stalk_ = list(self.s.pop(stalk))
                stalk_[0],l1.left,l2 = stalk_[0].common_substring(stalk)
                branch.add(l1.left,stalk_[1].copy())
                stalk_[1].update(reads)
                branch.add(l2,reads)
                stalk_ = tuple(stalk_)
                self[stalk_[0]] = branch
                self.s[stalk_[0]] = stalk_
            else:
                stalk_ = list(self.s.pop(stalk))
                branch = self.pop(stalk)
                stalk_[0],bstalk,stalk = stalk_[0].common_substring(stalk)
                if len(bstalk):
                    br = Branch()
                    br[bstalk] = branch 
                    br.s[bstalk] = (bstalk,stalk_[1].copy())
                    br.add(stalk,reads)
                    self[stalk_[0]] = br
                else: 
                    branch.add(stalk,reads)
                stalk_[1].update(reads)
                stalk_ = tuple(stalk_)
                if not len(bstalk): self[stalk_[0]] = branch
                self.s[stalk_[0]] = stalk_
        else:
            if type(stalk) == str: stalk = Stalk(stalk)
            self.s[stalk] = (stalk,reads)
            self[stalk] = Leaf(stalk)

'''
DESCRIPTION
    an object which constructs a suffix trie out of fragments of a sequence and can traverse 
    the trie to resconstruct some target sequence
INPUT
    reads | a list of strings which overlap and are fragments of a longer sequence
'''
class Sequitur:
    def __init__(self,reads,k_min=3,**kwargs):
        if "correct_sequence" in kwargs: self.correct_sequence = kwargs["correct_sequence"]
        self.branch = Branch()
        self.reads = reads#self.remove_containments(reads)
        self.k_min = k_min
        for read in reads:
            for i in range(len(read)):
                if len(read[i:]) < self.k_min: continue 
                self.branch.add(Stalk(read[i:]),{read})
        if "len" in kwargs or ("assemble" in kwargs and kwargs["assemble"]): self.phase1(**kwargs)

    def phase1(self,**kwargs): 
        if len(self.reads) == 1: 
            self.sequence = self.reads[0]
            return True
        else:
            if "k_min_add" not in kwargs: kwargs["k_min_add"] = 0   
            extensions = {}
            stalks = self.branch.b.keys()
            for stalk in stalks: self.longest_common_substring(self.branch,stalk,[stalk.stalk],extensions)
            k_max = max(extensions.keys())
            i = 0
            overlaps = {}
            if "len" in kwargs and kwargs["len"] == len(self.reads):
                if self.k_min + kwargs["k_min_add"] < k_max: kwargs["k_min_add"] += 1
                else:
                    if kwargs["assemble"]: kwargs["assemble"] = False
                    for v in extensions.values():
                        for k_, v_ in v.items():
                            for read in v_['endswith']:
                                if read not in overlaps: overlaps[read] = set()
                                for extension in v_['is_in']:
                                    for _ in range(extension.count(k_)):
                                        pre = extension[:extension.find(k_)+len(k_)]
                                        suf = extension[extension.find(k_)+len(k_):]
                                        if read.endswith(pre): overlaps[read].add((extension,suf))
                    if "biphasic" in kwargs and kwargs["biphasic"]:
                        if self.phase2(overlaps,**kwargs): return True
                        else: 
                            self.sequence = self.reads
                            return False
                    else:
                        self.sequence = self.reads
                        return True
            else:
                kwargs["k_min_add"] = 0 
                kwargs["len"] = len(self.reads)
            for read in self.reads:
                while min(k_max,len(read)-1) - i > self.k_min + kwargs["k_min_add"]:
                    if read[:min(k_max,len(read)-1)-i] in extensions[min(k_max,len(read)-1)-i]:
                        if len(extensions[min(k_max,len(read)-1)-i][read[:min(k_max,len(read)-1)-i]]['endswith']) > 1: i += 1
                        else: break
                    else: i+=1
                if (read[:min(k_max,len(read)-1)-i] in extensions[min(k_max,len(read)-1)-i] and len(extensions[min(k_max,len(read)-1)-i][read[:min(k_max,len(read)-1)-i]]['endswith']) > 1)\
                or read[:min(k_max,len(read)-1)-i] not in extensions[min(k_max,len(read)-1)-i] or read not in extensions[min(k_max,len(read)-1)-i][read[:min(k_max,len(read)-1)-i]]['is_in']:
                    i = 0
                    continue
                if list(extensions[min(k_max,len(read)-1)-i][read[:min(k_max,len(read)-1)-i]]['endswith'])[0] not in overlaps: overlaps[list(extensions[min(k_max,len(read)-1)-i][read[:min(k_max,len(read)-1)-i]]['endswith'])[0]] = (read[:min(k_max,len(read)-1)-i],read,read[min(k_max,len(read)-1)-i:])
                else: overlaps[list(extensions[min(k_max,len(read)-1)-i][read[:min(k_max,len(read)-1)-i]]['endswith'])[0]] = ('','','')
                i = 0
            if len(overlaps):
                overlaps = list(overlaps.items())
                overlaps.sort(key=lambda e: len(e[1][0]),reverse=True)
                overlaps = dict(overlaps)
                key = list(overlaps.keys())[0]
                seq = key
                self.reads.remove(key)
                while key in overlaps:
                    if len(overlaps[key][0]) < sum(len(o[0]) for o in overlaps.values())/len(overlaps): break
                    seq += overlaps[key][2]
                    key = overlaps[key][1]
                    if key not in self.reads or not len(key): break
                    self.reads.remove(key)
                self.reads += [seq]
            self.__init__(self.reads,self.k_min,**kwargs)
        
    def phase2(self,overlaps,**kwargs):
        for read,extensions in overlaps.items():
            if not len(extensions): continue
            self.reads.remove(read)
            for extension in extensions:
                self.reads.remove(extension[0])
                kwargs["len"] = len(self.reads)+1
                self.__init__(self.reads+[read+extension[1]],self.k_min,**kwargs)
                if self.phase1(**kwargs): return True
                else: 
                    self.reads.remove(read+extension[1])
                    self.reads += [extension[0]]
            self.reads += [read]
        return False
        
    def longest_common_substring(self,branch,stalk,substring,extensions):
        if branch.__is_shallow__():
            if len(branch.s[stalk][1]) > 1:
                if len(''.join(substring)) not in extensions: extensions[len(''.join(substring))] = {}
                if ''.join(substring) not in extensions[len(''.join(substring))]: extensions[len(''.join(substring))][''.join(substring)] = {'endswith':set(),'is_in':set()}
                for read in branch.s[stalk][1]:
                    if read.endswith(''.join(substring)): extensions[len(''.join(substring))][''.join(substring)]['endswith'].add(read)
                    else: extensions[len(''.join(substring))][''.join(substring)]['is_in'].add(read)
            return extensions
        if type(branch.__traverse__(stalk.stalk)) is Leaf:
            if len(branch.s[stalk][1]) > 1:
                if len(''.join(substring)) not in extensions: extensions[len(''.join(substring))] = {}
                if ''.join(substring) not in extensions[len(''.join(substring))]: extensions[len(''.join(substring))][''.join(substring)] = {'endswith':set(),'is_in':set()}
                for read in branch.s[stalk][1]:
                    if read.endswith(''.join(substring)): extensions[len(''.join(substring))][''.join(substring)]['endswith'].add(read)
                    else: extensions[len(''.join(substring))][''.join(substring)]['is_in'].add(read)
            return extensions
        for c in branch.__traverse__(stalk.stalk).b:
            if c.stalk != '$': extensions = self.longest_common_substring(branch.__traverse__(stalk.stalk),c,substring+[c.stalk],extensions)
            else: 
                if len(branch.s[stalk][1]) > 1:
                    if len(''.join(substring)) not in extensions: extensions[len(''.join(substring))] = {}
                    if ''.join(substring) not in extensions[len(''.join(substring))]: extensions[len(''.join(substring))][''.join(substring)] = {'endswith':set(),'is_in':set()}
                    for read in branch.s[stalk][1]:
                        if read.endswith(''.join(substring)): extensions[len(''.join(substring))][''.join(substring)]['endswith'].add(read)
                        else: extensions[len(''.join(substring))][''.join(substring)]['is_in'].add(read)
        return extensions
    
    def remove_containments(self,reads):
        i = 0
        r = set()
        m = max([len(r) for r in reads])
        b = False
        while i < len(reads):
            if len(reads[i]) == m: 
                r.add(reads[i])
                i+=1
                continue
            for r_ in r:
                if reads[i] in r_: 
                    b = True
                    break
            if b: 
                b = False
                i+=1
                continue
            r.add(reads[i])
            i+=1
        return list(r)
