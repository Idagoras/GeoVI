import json


class CheckInData:
    def __init__(self, user_id, time_stamp, latitude, longitude, loc_id):
        self.user_id = user_id
        self.time_stamp = time_stamp
        self.latitude = latitude
        self.longitude = longitude
        self.loc_id = loc_id


def is_float(value):
    try:
        float_value = float(value)
        return True
    except ValueError:
        return False


with open('/Users/wuwei/workspace/GeoVI/dataset/loc-brightkite_totalCheckins.txt', 'r') as file:
    # 逐行读取文件内容
    for line in file:
        words = line.split()
        try:
            m_user_id = int(words[0])
            m_time_stamp = words[1]
            m_latitude = float(words[2])
            m_longitude = float(words[3])
            m_loc_id = words[4]
            if 48.9086 >= m_latitude >= 48.8091:
                if 2.4688 >= m_longitude >= 2.225:
                    with open('Paris.txt', 'a') as paris_checkin_file:
                        paris_checkin_file.write(line)
                        print("write to paris.txt\n")
            if -74.2557 <= m_longitude <= -73.7682:
                if 40.5836 <= m_latitude <= 40.8132:
                    with open('New_York.txt', 'a') as new_york_checkin_file:
                        new_york_checkin_file.write(line)
                        print("write to new york.txt\n")
        except ValueError:
            print("can not parse the string line\n")
        except IndexError:
            print("format error\n")
