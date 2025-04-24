# color
cell_type_colors = {
   'Str': '#FEE9B4', 'NI': '#ED5A65', 
    'Meso': '#EF859B', 'CHGB+ CALB2+ nTh': '#9F5F9F', 
    'Th': '#C2C1E0', 'WM': '#B9E8EA', 'NIVL': '#FF7F00',
    'SGF': '#98df8a', 'ZIC1+ ZIC2+ nTh': '#F3BAA5', 
    'SGC': '#68945c', 'Vasc': '#C6808C', 'SOp': '#9bb496', 
    'HSPB7+ nTh': '#EAADEA', 'Arco': '#FFC0CB',
    'Hyper': '#FFCCCC', 'Ento': '#cde6c7',
"CGL": "#F5F7D0",
"CML": "#ffd92f",
"NF": "#f781bf",
"BS": "#8da0cb",
"CPL": "#ffec61",
"Hp": "#C0D9D9",
"TeO": "#9370DB", 'NC': '#F08080', 'Field L': '#D9D9F3', 'SGP': '#67ADB7', 'Imc/Ipc': '#D3E2B7',
    'Field L': '#D9D9F3',"NCC": "#ff461f",
    "NCM": "#f47983", 
    "NCL": "#F35336" 
}

#Plot
spot_size = 150
sc.pl.spatial(adata, color='simplified_annotated_spatial_leiden', spot_size=spot_size, palette=cell_type_colors)