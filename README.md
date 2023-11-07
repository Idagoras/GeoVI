# GeoVI  ![Static Badge](https://img.shields.io/badge/build-passed-orange)

Geographic indistinguishability mechanism under Voronoi diagram partitioning

---

使用维诺图划分的地理不可区分性机制

## 主要内容

+ ✅读取 OpenStreetMap 通用格式 `.osm ` 的地图文件，并解析成无向/有向图。
+ ✅选择地图中的感兴趣点集合，根据感兴趣点集合生成相应的维诺图（默认距离度量是欧几里得距离）

  + ✅经纬度坐标点（ WGS84 ）到 UTM 投影坐标系下的二维坐标点转换
  + ✅维诺图的可视化
+ ✅计算维诺图中 cell 之间的最短距离
+ ⭕为每个 cell 的添加语义敏感度

  + ⭕语义敏感度的计算方法
  + ⭕语义敏感度的正确性
+ ⭕合并 cell 得到满足隐藏区域要求的语义敏感值的隐藏区域

  + ⭕合并的背包问题算法
+ ⭕为每个隐藏区域的语义敏感度添加差分隐私噪声

  + ☑拉普拉斯噪声
+ ⭕实现地点扰动的地理图指数机制

  + ✅地理图指数机制
  + ⭕扰动区域内 K - 匿名机制的实现
+ ⭕服务质量的优化算法实现

## 相关依赖库

+ Boost
+ Osmium
+ Proj
+ OpenSceneGraph
+ ....
