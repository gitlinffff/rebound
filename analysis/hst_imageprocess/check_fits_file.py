from astropy.io import fits
from datetime import datetime, timedelta

"""
HST
"""
# Load FITS image
with fits.open('/home/linfel/linfel_data/hst_raw_JianyangLi/17289/stack_all_long.fits') as hdul:
    hst_data = hdul[0].data  # assuming image is in HDU 0
    header = hdul[0].header

# Extract metadata
orientat = header.get('ORIENTAT', 'N/A')
utc_mid = header.get('UTC-MID', 'N/A')

# the time of the hst image
hst_time = datetime.strptime(utc_mid, "%Y-%m-%dT%H:%M:%S.%f")

# Time of the impact
impact_datetime = "2022-09-26T23:14:24.1830"
impact_time = datetime.strptime(impact_datetime, "%Y-%m-%dT%H:%M:%S.%f")

# calculate the time difference
time_difference = hst_time - impact_time
days_after_impact = time_difference.total_seconds() / 86400.
print(f"Time of HST data: T0+{days_after_impact:.2f} days")
