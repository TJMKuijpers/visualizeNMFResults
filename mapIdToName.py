import pandas as pd


def mapidtoname(idstomap=None, database=None, idcolumn=None):
    idsMappedToDataBase=database[database[idcolumn].isin(idstomap)]
    return idsMappedToDataBase


def mapidtocpgvariable(idstomap=None, database=None, idcolumn=None, columntoreturn=None):
    # get the gene ids that are mapped to the database
    idstoextract = database[database[idcolumn].isin(idstomap)]
    # return the CpG variable for that gene
    idstoreturn = idstoextract[columntoreturn]
    return idstoreturn


def setthedatabase(pathtodatabase=None, seperator=None):
    database = pd.read_csv(pathtodatabase, sep=seperator)
    return database


def mapgenenamestoctd(genenames=None, ctddatabase=None, columnid=None):
    idsmappedtoctd = ctddatabase[ctddatabase[columnid].isin(genenames)]
    return idsmappedtoctd
