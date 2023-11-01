
class fastaSeqFile:
    def __init__(self, fileName):
        try:
            self.file = file(fileName)
        except:
            print 'cannot open "%s"' % fileName
            raise
        self.fileName = fileName
        self.iterMode = 0

        # skip leading junk
        self.fileIt = iter(self.file)
        self.defLine = ''
        for l in self.fileIt:
            if '>' == l[0]:
                self.defLine = l
                break

    def __iter__(self):
        # new iterator, new object.
        if self.iterMode: return self
        else:
            newSelf = fastaSeqFile(self.fileName)
            newSelf.iterMode = 1
            return newSelf

    def next(self):
        if '' == self.defLine: raise StopIteration
        s = [self.defLine]
        self.defLine = ''
        for l in self.fileIt:
            if '>' == l[0]:
                self.defLine = l
                break
            s.append(l)
        return s
