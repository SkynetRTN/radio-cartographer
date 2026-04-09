from astropy.io import fits
import sys

def check_timestamps(file_path):
    print(f"Opening {file_path}...")
    try:
        hdul = fits.open(file_path)
        data = hdul[1].data
        header = hdul[1].header
        # 'DATE-OBS' contains the observation timestamps
        #print(header)
        time = data['MJD']
        #print(time[:200])
        total_rows = len(time)
        dupes = sum(time[i] == time[i-1] for i in range(1, total_rows))
        
        print(f"Total rows: {total_rows}")
        print(f"Duplicate consecutive timestamps: {dupes}")
        
    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    target_file = '/skynet/radio-cartographer/testing/test_files/0144845.fits'
    if len(sys.argv) > 1:
        target_file = sys.argv[1]
        
    check_timestamps(target_file)
