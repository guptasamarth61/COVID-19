import pandas as pd
df = pd.read_csv("files\Italy_Change in Mobility.csv",header = 0 ,names=['c_id','country','sr1','sr2','date','retail','grocery','parks','transit','workplaces','residential'])
df = df.groupby('date').mean()
df.to_csv("Italy Mobility.csv")
