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
    def OpenSection(this,name,num):
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
    def CloseSection():
        return IOSection.CloseSection(this.handle)
    def ReadVar(this,name):
        return IOSection.ReadVar(this.handle,name)





a=IOSectionClass()
a.OpenFile("grr.h5")
a.OpenSection("Energies")
print len(a.ReadVar("PotentialEnergy"))
