#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>

#include <xoos/types/float.h>
#include <xoos/types/int.h>

#include "fastcluster.h"

namespace xoos::cnc::segmentation {

/**
 * @brief a point in 2D space, can be extended to n-dimensions
 * A data structure for calculating distance between data points
 */
struct Point {
  f64 x;
  f64 y;

  Point(f64 xx, f64 yy) : x(xx), y(yy) {
  }
};

/**
 * @brief Compute the Euclidean distance between two points
 */
static inline f64 GetDistance(const Point& p1, const Point& p2) {
  f64 distance = std::sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
  return distance;
}

/**
 * @brief Wrapper function of the hclust_fast
 * @return Vector of cluster id from 0, 1, ..., ncluster - 1
 */
std::vector<s32> Hclust(const std::vector<Point>& points, s32 hclust_method, const f64 h) {
  auto npoints = static_cast<s32>(points.size());

  std::vector<f64> distmat;
  std::vector<s32> merge;
  std::vector<f64> height;
  std::vector<s32> labels;

  distmat.resize(npoints * (npoints - 1) / 2);
  merge.resize(2 * (npoints - 1));
  height.resize(npoints - 1);
  labels.resize(npoints);

  // Compute the upper triangle of the distance matrix excluding diagnol values
  s32 k = 0;
  for (s32 i = 0; i < npoints; i++) {
    for (s32 j = i + 1; j < npoints; j++) {
      distmat[k++] = GetDistance(points[i], points[j]);
    }
  }

  if (npoints > 1) {
    hclust_fast(npoints, distmat.data(), hclust_method, merge.data(), height.data());
    cutree_cdist(npoints, merge.data(), height.data(), h, labels.data());
    return labels;
  } else if (npoints == 1) {
    return {0};
  } else {
    throw std::runtime_error("hclust: npoints is 0!");
  }
}

}  // namespace xoos::cnc::segmentation
