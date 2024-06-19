import pandas as pd


def getSamples(iptFile):
    df = pd.read_table(iptFile)
    return df.SampleName.tolist()


def getPair(iptFile):
    df = pd.read_table(iptFile)
    dct = {}
    for row in df.to_dict(orient="records"):
        if isinstance(row["ControlName"], str):
            dct[row["SampleName"]] = [row["SampleName"], row["ControlName"]]
    return dct
