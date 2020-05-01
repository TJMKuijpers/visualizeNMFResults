import pandas as pd


def mapidtoname(idstomap=None, database=None,idColumn=None):
    idsMappedToDataBase=database[database[idColumn].isin(idstomap)]
    return idsMappedToDataBase


def setthedatabase(pathtodatabase=None, seperator=None):
    database = pd.read_csv(pathtodatabase, sep=seperator)
    return database


def mapgenenamestoctd(genenames=None,ctdDataBase=None,columnId=None):
    idsMappedToCtd=ctdDataBase[ctdDataBase[columnId].isin(genenames)]
    return idsMappedToCtd
