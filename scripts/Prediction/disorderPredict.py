#!/usr/bin/env python3

import urllib.request
import json
import pandas as pd 


def getDisorder(protein= 'A4V4D2'):
    
    try:  
        # Define request
        acceptHeader = 'application/json' # text/csv and text/plain supported
        # request = urllib.request.Request("http://mobidb.bio.unipd.it/ws/P04050/disorder", headers={"Accept" : acceptHeader})
    
        path = "http://mobidb.bio.unipd.it/ws/" + protein + "/consensus"
        request = urllib.request.Request(path, headers={"Accept" : acceptHeader})
    
        # Send request
        response = urllib.request.urlopen(request)
    
        # Parse JSON response di Python dict
        data = json.load(response)

        data["mobidb_consensus"]['disorder']['predictors'][1]['method']
        mobidb_lite = data["mobidb_consensus"]['disorder']['predictors'][1]
        return mobidb_lite['regions']
    except (IndexError, TypeError):
        return None
    
def formatDisorder(l):
    if l == []:
        return ''
    if l:
        flattened_list =  ",".join([":".join([str(l) for l in sublist[0:2]]) for sublist in l])
        return flattened_list
    return l 
    
df = pd.read_csv("../../output/adjusted/fly/mergedPolyAaDf.csv")
df = df.loc[1:15]
df['disorderedPred'] = df.apply(lambda row : formatDisorder(getDisorder(row.uniprotsptrembl)), axis=1)

df.to_csv("../../output/adjusted/fly/mergedPolyAaDfDisorder.csv")
