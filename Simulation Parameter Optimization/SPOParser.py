from dataclasses import Field
import re
import SPO
# TODO:
# [x] Config File Parser
# [ ] Ensemble Writer
# [ ] Log File Parser



#Parser should be simple, put complexity in the file spec
class SPOEnsembleLogParser:
    def __init__(self,fileName):
        self.file = fileinput.input(fileName,inplace=True)

    def __del__(self):
        self.file.close()
    
    def parseEnsembleLog(self,ensembleID):
        finished = True
        for line in self.file:
            #look for lines that end in "Running " 
            runningJob = re.match("Run (\d+) of \d+ Running ",line)
            if runningJob:
                #we found a job that is still running
                if int(runningJob.group(1)) == ensembleID:
                    # we are the responsible of that job,
                    # update that is is finished
                    line = line.replace("Running","Finished")
                else:
                    finished = False
            print(line,end='')
        return finished

#the parser works by filling out a file spec
# a file Spec is made out of fieldSpecs.
# field specs either are header and a rule to parse the lines after
# or 

# two options here, inherite from fileSpec every time I need a thing
# or make a 
class FileSpecFactory:
    def __init__(self):
        self.floatMatch = "([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
        
    def logFileSpec(self):
        logField = FieldSpec(format="Step\s+(\d)\s+{(.*)}\s+Residue:"+self.floatMatch)
        lastLogField = FieldSpec(format="Step\s+(\d)\s+{(.*)}\s+")
        logSpec = FileSpec()
        return logSpec
        
    def configFileSpec(self):
        #create the empty spec
        configSpec = HeaderFieldFileSpec()
        # the name field take everything on the next line 
        # so we can leave it blank
        nameField = HeaderFieldSpec("Name:")
        configSpec.addFieldSpec(nameField)

        simulationField = HeaderFieldSpec("Simulation:")
        configSpec.addFieldSpec(simulationField)

        partitionField = HeaderFieldSpec("Partition:")
        configSpec.addFieldSpec(partitionField)
        
        extraScriptCommandsField = HeaderFieldSpec("Extra Script Commands:",multiLine = True)
        configSpec.addFieldSpec(extraScriptCommandsField)

        methodField = HeaderFieldSpec("Method:")
        configSpec.addFieldSpec(methodField)
        # dataFileField matches file names
        # look for some not white space then a period then some word characters
        fileNameRegex = "(\S*\.w+)" 
        dataFileField = HeaderFieldSpec("Data:",\
            format=fileNameRegex + "\s+" + fileNameRegex,multiLine=True)
        configSpec.addFieldSpec(dataFileField)
        # parameters look for alphanumerics and underscores followed by a colon
        # then a float
        paramField = HeaderFieldSpec("Parameters:",\
            format = "(.*):\s*"+self.floatMatch,multiLine = True)
        configSpec.addFieldSpec(paramField)

        # now the more complicated fields
        # runs on need to make sure we have Either Desktop or HPCC and might 
        # have a number after it. So we define some post processing to
        # check these things
        def runsOnPostProcess(self):
            # check if we match "Desktop"
            runsOnMatch = re.match("Desktop",self.value)
            if runsOnMatch:
                self.value = [SPO.SPORunsOn.Desktop]
            else:
            # check if we match HPCC
                line = self.value
                runsOnMatch = re.match("HPCC",line)
                if runsOnMatch:
                    self.value = [SPO.SPORunsOn.HPCC]
                    #check if we have a number after
                    runsOnSplit = line.split(' ')
                    if len(runsOnSplit)>1:
                        # if so append it
                        self.value.append(int(runsOnSplit[1]))
                    else:
                        # if not append 1
                        self.value.append(1)
                else:
                    raise Exception("Field 'Runs On' is formated incorrectly.\
                        Expected 'Desktop' or 'HPCC' possibly followed by a number")
        
        runsOnField = HeaderFieldSpec("Runs On:",postProcessing=runsOnPostProcess)
        configSpec.addFieldSpec(runsOnField)

        return configSpec

class FileSpec:
    # a Plane old file that tries to read a single field
    def __init__(self):
        self.fields = []
    def addFieldSpec(self,field):
        self.fields.append(field)
    def readFile(self,fileName):
        #for each line try each field
        parser = SPOFileParser(fileName)
        # we'll rely on the parser raising a stop iteration when its done
        # with the file
        while parser:
            self.findAndFillSpec(parser)

        return self.makeValues()

    def findAndFillSpec(self,parser):
        # loop though all the fields until we get one that can read the line
        for field in self.fields:
            if field.read(parser):
                break

    def makeValues(self):
        return list(field.value for field in self.fields)
    

class HeaderFieldFileSpec(FileSpec):
    def __init__(self):
        self.fields = []
    def addFieldSpec(self,field):
        self.fields[field.header]=field

    def findAndFillSpec(self,parser):
        line = next(parser)
        # Find the field with a header that matches and read it
        if line in self.fields.keys():
            if not self.fields[line].read(parser):
                raise Exception("Field %s is formatted incorrectly could not\
                                match regex %s"%(self.fields[line].header,self.fields[line].format))
        else:
            raise Exception("Unrecognized header in file %s"%parser.fileName)

    def makeValues(self):
        #now return a dictionary of file headers and their values
        # for each field in the list of fields, make a dictionary of {header:value}
        return dict(list((field.header,field.value) for field in self.fields))


class FieldSpec:
    def __init__(self,format = "(.*)",multiLine = True,postProcessing=None):
        self.format = format
        self.multiLine = multiLine
        self.postProcessing = None

    def matchLine(self,parser):
        # assume what we've been given is a regular expression
        # find the matches and report all of them in a list
        line = next(parser)
        matches = re.match(format,line)
        if matches:
            # exclude the first element of match as 
            # it contains the whole match
            return matches[1::]
        return None


    def read(self, parser):
        # Tries to match the current line to the pattern return true if we have 
        # a match and false if we don't
        # check if our format is a callable if so call it
        if callable(self.format):
            self.format(parser)
        elif isinstance(self.format,str):
            # always read the first line
            self.value = self.matchLine(parser)
            # make sure we could match the line
            if not self.value:
                return False

            if self.multiLine:
                # for multi line field value should be 
                # a list of value on the lines

                self.value = [self.value] 
                nextVal = self.matchLine(parser)
                while nextVal:
                    self.value.append(nextVal)
                    nextVal = self.matchLine(parser)
                # This kind of read will always go one past the end so we need to rewind
                parser.rewind()
                
        # call postprocessing if we need to
        if self.postProcessing is not None:
            self.postProcessing()
        return True
        

class HeaderFieldSpec():
    def __init__(self,header,format = "(.*)",postProcessing=None,multiLine = False):
        self.header = header
        self.format = format #default value for this takes the entire line
        self.value = None
        self.postProcessing = postProcessing
        self.multiLine = multiLine



class SPOFileParser:
    def __init__(self,fileName, mode = None):
        self.file = open(fileName,"r+")
        self.fileName = fileName
        #keep a list of line positions so we can rewind
        self.lastLinePos = []
        self.lastLinePos.append(self.file.tell())
    def __del__(self):
        self.file.close()

    def parseLogFile(self):
        step = 0
        params = []
        residues = []
        #skip header
        next(self)
        next(self)
        next(self)
        for line in self:
            logData = re.match("Step\s+(\d)\s+{(.*)}.*",line)
            if logData:
                step = int(logData.group(1))
                paramList = []
                paramPairs = logData.group(2).split(",")
                for param in paramPairs:
                    (name,value) = param.split(':')
                    paramList.append([name,float(value)])
                params.append(paramList)
                # see if there is a residue
                resMatch = re.match(".*Residue:"+self.floatMatch,line)
                if resMatch:
                    residues.append(float(resMatch.group(1)))

            else:
                raise Exception("Error reading log file")
        return(params,residues,step)

    def __next__(self):
        #grabs the next line that isn't just a new line or comment ('#')
        lastLine = self.file.tell()
        self.lastLinePos.append(lastLine)
        
        line = ""
        while len(line) == 0:
            line = self.file.readline()
            if self.file.tell() == lastLine:
                raise StopIteration
            line = self._trimLine(line)
            lastLine = self.file.tell()

        return line
        
    def __iter__(self):
        return self
    def _trimLine(self, line):
        #look for the comment symbol
        comment = re.search("(.*)#",line)
        # if its found
        if comment:
            line = comment.group(1) # grab the first group as the line
        line = line.strip() #remove any leading or trailing whitespace
        return line
    def rewind(self):
        self.file.seek(self.lastLinePos.pop())
