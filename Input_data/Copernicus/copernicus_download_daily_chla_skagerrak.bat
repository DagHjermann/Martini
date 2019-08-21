:: Downloads year-by-year data from
:: North Atlantic Chlorophyll (Copernicus-GlobColour) from Satellite Observations: Daily Interpolated (Reprocessed from 1997) 
:: OCEANCOLOUR_ATL_CHL_L4_REP_OBSERVATIONS_009_098
:: Data are downloaded to folder C:/Data/temp
:: How did I make this file? 
::   First, install motu client by python -m pip install motuclient (see http://marine.copernicus.eu/faq/what-are-the-motu-and-python-requirements/ )
::   Went to the data set page (http://marine.copernicus.eu/services-portfolio/access-to-products/?option=com_csw&view=details&product_id=OCEANCOLOUR_ATL_CHL_L4_REP_OBSERVATIONS_009_098) 
::   1) Click download, 2)  pick area and times, 3) CLick "download options", 4) Click "View script"
::
:: Looping in batch files: https://stackoverflow.com/a/2591782/1734247 

for /l %%x in (2001,1,2018) do (
  python -m motuclient --motu http://my.cmems-du.eu/motu-web/Motu --service-id OCEANCOLOUR_ATL_CHL_L4_REP_OBSERVATIONS_009_098-TDS --product-id dataset-oc-atl-chl-multi-l4-oi_1km_daily-rep-v02 --longitude-min 6.75 --longitude-max 12.994793891906738 --latitude-min 57.5 --latitude-max 60 --date-min "%%x-01-01 00:00:00" --date-max "%%x-12-31 00:00:00" --variable CHL --variable CHL_error --out-dir "C:/Data/temp" --out-name "dataset-oc-atl-chl-multi-l4-oi_1km_daily-rep-v02_%%x.nc" --user Dhjermann --pwd "decticus22&"
  )
