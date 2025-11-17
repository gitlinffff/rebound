from datetime import datetime

# Define the two time strings
time_str_1 = "2022-09-26 23:14:24.1830"
time_str_2 = "2022-11-30 09:49:36.3960"

# Define the format string to correctly parse the times, including milliseconds/microseconds
time_format = "%Y-%m-%d %H:%M:%S.%f"

# Convert the strings to datetime objects
start_time = datetime.strptime(time_str_1, time_format)
end_time = datetime.strptime(time_str_2, time_format)

# Calculate the time difference (timedelta)
time_difference = end_time - start_time

# Convert the time difference to the total number of seconds
total_seconds = time_difference.total_seconds()

# Print the result
print(f"Start Time: {start_time}")
print(f"End Time:   {end_time}")
print(f"Time Difference (timedelta): {time_difference}")
print(f"Time Difference in **Seconds**: {total_seconds}")
print(f"Time Difference in **Days**: {total_seconds/86400.}")
