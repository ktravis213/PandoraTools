# PandoraTools

GetPandorafromRSIG.py uses calls to EPA's RSIG tool to get Pandora sky-scan and direct-sun NO2 for a given date-range across CONUS

MatchPandoraandTEMPO.py reads TEMPO data and matches it to Pandora data within 5km and 30min. It assumes an already processed Pandora dataframe from GetPandorafromRSIG.py
