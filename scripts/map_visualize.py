import folium
import xml.etree.ElementTree as ET
import random


def generate_random_colors(k):
    colors = []
    for _ in range(k):
        # 生成随机的六位十六进制颜色值
        color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
        colors.append(color)
    return colors


tree = ET.parse('/Users/wuwei/workspace/GeoVI/build/Paris_test.xml')
root = tree.getroot()

colors = ['lightgreen', 'green', 'darkpurple', 'orange', 'white', 'gray', 'pink', 'lightblue', 'darkgreen', 'cadetblue',
          'beige', 'darkred', 'lightred', 'lightgray', 'blue', 'purple', 'black', 'red', 'darkblue']

index = 0
k = 0
m = folium.Map(location=[(48.8531000 + 48.8439000) / 2, (2.346300 + 2.3317000) / 2], zoom_start=18)
clusters = root.findall('.//cluster')
cluster_size = len(clusters)
random_colors = generate_random_colors(cluster_size)
for cluster in clusters:
    for location in cluster.findall('.//location'):
        latitude = float(location.find('latitude').text)
        longitude = float(location.find('longitude').text)
        print(f"cluster {index} location latitude = {latitude} longitude = {longitude}")
        folium.Marker(location=[latitude, longitude], popup=f"cluster{index}",
                      icon=folium.Icon(color=colors[k % len(colors)])).add_to(m)
    index += 1
    k += 1

m.save('map.html')
