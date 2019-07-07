#!/usr/bin/env python3

import urllib.request
import json
import pandas

def getDisorder(protein= 'A4V4D2'):

    # Define request
    acceptHeader = 'application/json' # text/csv and text/plain supported
    # request = urllib.request.Request("http://mobidb.bio.unipd.it/ws/P04050/disorder", headers={"Accept" : acceptHeader})

    path = "http://mobidb.bio.unipd.it/ws/" + protein + "/consensus"
    request = urllib.request.Request(path, headers={"Accept" : acceptHeader})


    # Send request
    response = urllib.request.urlopen(request)

    # Parse JSON response di Python dict
    data = json.load(response)

    try:  
        data["mobidb_consensus"]['disorder']['predictors'][1]['method']
        mobidb_lite = data["mobidb_consensus"]['disorder']['predictors'][1]
        return mobidb_lite['regions']
    except IndexError:
        return [None]
    
# Figure
# handle data
# =============================================================================
# print(data)
# print(disorder)
# print(db)
# =============================================================================


#print(getDisorder("Q9W5X6")) # bad


# bad X1
#protein = "Q9W5W6" # bad
## Define request
#acceptHeader = 'application/json' # text/csv and text/plain supported
## request = urllib.request.Request("http://mobidb.bio.unipd.it/ws/P04050/disorder", headers={"Accept" : acceptHeader})
#
#path = "http://mobidb.bio.unipd.it/ws/" + protein + "/consensus"
#request = urllib.request.Request(path, headers={"Accept" : acceptHeader})
#
#
## Send request
#response = urllib.request.urlopen(request)
#
## Parse JSON response di Python dict
#data = json.load(response)
#
#data["mobidb_consensus"]['disorder']['predictors'][1]['method']
#mobidb_lite = data["mobidb_consensus"]['disorder']['predictors'][1]
#print(mobidb_lite)
# Figure