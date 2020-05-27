def getchromsomelocations(cpgsites=None,dnainfo=None,columnname=None):
    cpgInformation={}
    for x in cpgsites:
        # map the CpG sites against the Sites in the database dnainfo
        tmp_info=dnainfo[dnainfo[columnname].isin(x)]
        cpgInformation['cluster'+str(cpgsites.index(x))]=tmp_info.loc[:,['TargetID','CHR']]
    return cpgInformation

