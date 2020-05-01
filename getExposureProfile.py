def getexposureprofile(inputData=None,idSampleColumn=None,clusterData=None,columnsToSelect=None):
    popProfile=inputData[inputData[idSampleColumn].isin(clusterData)]
    if columnsToSelect==None:
        popProfile = popProfile
    else:
        popProfile=popProfile[:,columnsToSelect]
    return popProfile
