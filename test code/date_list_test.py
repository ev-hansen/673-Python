from datetime import datetime, timedelta


def date_range(start, end):
    delta = end - start  # as timedelta
    days = [start + timedelta(days=i) for i in range(delta.days + 1)]
    return days


given_time_s_str = "2013-03-17/00:00:00"
given_time_e_str = "2013-05-18/00:00:00"

time_s_datetime = datetime.strptime(given_time_s_str, '%Y-%m-%d/%H:%M:%S')
time_e_datetime = datetime.strptime(given_time_e_str, '%Y-%m-%d/%H:%M:%S')

print(date_range(time_s_datetime, time_e_datetime))
