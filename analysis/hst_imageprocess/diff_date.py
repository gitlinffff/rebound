from datetime import datetime

def time_diff(time_str_1, time_str_2):
	# Define the format string to correctly parse the times, including milliseconds/microseconds
	time_format = "%Y-%m-%dT%H:%M:%S.%f"

	# Convert the strings to datetime objects
	start_time = datetime.strptime(time_str_1, time_format)
	end_time = datetime.strptime(time_str_2, time_format)

	# Calculate the time difference (timedelta)
	time_difference = end_time - start_time

	# Convert the time difference to the total number of seconds
	total_seconds = time_difference.total_seconds()

	# Print the result
	print(f"\n{'='*70}")
	print(f"Start Time: {start_time}")
	print(f"End Time:   {end_time}")
	print(f"Time Difference (timedelta): {time_difference}")
	print(f"Time Difference in **Seconds**: {total_seconds}")
	print(f"Time Difference in **Days**: {total_seconds/86400.}")


if __name__ == "__main__":
	t_str1 =  "2022-09-26T23:14:24.1830"
	t_str2 = ["2022-09-27T01:06:21.979",
						"2022-09-27T02:45:33.987",
						"2022-09-27T04:16:50.987",
						"2022-09-27T05:52:06.986",
						"2022-09-27T07:27:21.979",
						"2022-09-27T17:02:00.979",
						"2022-09-28T02:42:25.227",
						"2022-09-28T17:06:51.739",
						"2022-09-29T02:47:45.483",
						"2022-09-30T16:31:02.060",
						"2022-10-01T16:28:21.068",
						"2022-10-08T19:52:10.556",
						"2022-10-11T20:59:02.812",
						"2022-10-15T10:39:55.068",
						"2022-11-30T09:49:36.396",
						"2022-12-14T14:58:03.229",
						"2022-12-28T15:09:16.926"
						]
	for t2 in t_str2:
		time_diff(t_str1, t2)
