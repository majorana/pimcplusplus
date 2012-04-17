from HTMLgen import *

def NewTable():
     myTable=Table()
     myTable.border=1
     myTable.width='40%'
     myTable.column1_align='center'
     myTable.cell_align='center'
     return myTable

def BuildTable(rows,cols):
    myTable=NewTable()
    row=[]
    for colNum in range(0,cols):
        row.append([])
    myTable.body=[row]
    for rowNum in range(1,rows):
        row=[]
        for colNum in range(0,cols):
            row.append([])
        myTable.body.append(row)
    return myTable

