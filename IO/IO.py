import IOSection

class IOSectionClass:
    handle=0
    def __init__(this):
        this.handle=IOSection.New()
    def OpenFile(this,fileName):
        return IOSection.OpenFile(this.handle,fileName)
    def GetName(this):
        return IOSection.GetName(this.handle)
    def NewFile(this,fileName):
        return IOSection.NewFile(this.handle,fileName)
    def CloseFile(this):
        return IOSection.CloseFile(this.handle)
    def FlushFile(this):
        return IOSection.FlushFile(this.handle)
    def OpenSection2(this,name,num):
        return IOSection.OpenSectionNameNum(this.handle,name,num)
    def OpenSection(this,namenum):
        if type(namenum)==int:
            return IOSection.OpenSectionNum(this.handle,namenum)
        else:
            return IOSection.OpenSectionName(this.handle,namenum)
    def IncludeSection(this,name,fileName):
        return IOSection.IncludeSection(this.handle,name,fileName)
    def NewSection(this,name):
        return IOSection.IncludeSection(this.handle,name,fileName)
    def CloseSection(this):
        IOSection.CloseSection(this.handle)
    def ReadVar(this,name):
        return IOSection.ReadVar(this.handle,name)
    def CountSections(this):
        return IOSection.CountSections(this.handle)
    def CountSections2(this,name):
        return IOSection.CountSectionsName(this.handle,name)
    def CountVars(this):
        return IOSection.CountVars(this.handle)
    def GetVarName(this,num):
        return IOSection.GetVarName(this.handle,num)
    



#a=IOSectionClass()
#a.OpenFile("grr.h5")
#a.OpenSection("Energies")
#print len(a.ReadVar("PotentialEnergy"))
#a.CloseSection()
#a.OpenSection("ppPC");
#a.OpenSection("grid");
#r = a.ReadVar ("Points");
#print r

## b = IOSectionClass()
## b.OpenFile ("junk.dat")
## print b.ReadVar ("string_data")
## print b.ReadVar ("int_2");
## print b.ReadVar ("double_2");
## print b.ReadVar ("bool_2");
## print b.ReadVar ("string_2");

## print b.ReadVar ("int_3");
## print b.ReadVar ("double_3");
## print b.ReadVar ("bool_3");
## print b.ReadVar ("string_3");

## print 
## print b.ReadVar ("int_4");
## print 
## print b.ReadVar ("double_4");
## print 
## print b.ReadVar ("bool_4");
## print 
## print b.ReadVar ("string_4");
