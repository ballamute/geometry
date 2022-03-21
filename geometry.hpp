#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace Geometry {

class Vector;

class IShape;
class Point;
class Segment;
class Ray;
class Line;
class Circle;
class Polygon;

class Vector {
 public:
  using CoordinateType = int;
  using OperationResultType = long long;
  using MultiplierType = int;

  Vector();
  Vector(const CoordinateType& p_x, const CoordinateType& p_y);
  Vector(const Point& first_point, const Point& second_point);
  Vector(const Vector& vector);
  ~Vector() = default;

  Vector operator+=(const Vector& vector);

  Vector operator-=(const Vector& vector);

  Vector operator*=(const MultiplierType& mul);

  Vector operator-();

  Vector& operator=(const Vector& vector);

  OperationResultType GetScalarSquare() const;

  CoordinateType x;
  CoordinateType y;
};

class IShape {
 public:
  virtual ~IShape() = default;

  virtual IShape& Move(const Vector& vector) = 0;

  virtual bool ContainsPoint(const Point& point) const = 0;

  virtual bool CrossesSegment(const Segment& segment) const = 0;

  virtual IShape* Clone() const = 0;

  virtual std::string ToString() = 0;
};

class Point : public IShape {
 public:
  Point() = default;
  explicit Point(const Vector& vector);
  Point(const Vector::CoordinateType& p_x, const Vector::CoordinateType& p_y);

  IShape& Move(const Vector& vector) override;

  bool ContainsPoint(const Point& point) const override;

  bool CrossesSegment(const Segment& segment) const override;

  IShape* Clone() const override;

  std::string ToString() override;

  Point& operator=(const Point& point);

  Vector point_vector;
};

class Segment : public IShape {
 public:
  Segment() = default;
  Segment(const Point& first_point, const Point& second_point);
  Segment(const Segment& segment);

  IShape& Move(const Vector& vector) override;

  bool ContainsPoint(const Point& point) const override;

  bool CrossesSegment(const Segment& segment) const override;

  IShape* Clone() const override;

  std::string ToString() override;

  std::pair<Point, Point> GetPoints() const;
  Vector GetDirectionVector() const;

 private:
  Point first_point_;
  Point second_point_;
};

class Ray : public IShape {
 public:
  Ray() = default;
  Ray(const Ray& ray);
  Ray(const Point& start_point, const Vector& direction_vector);
  Ray(const Point& first_point, const Point& second_point);

  IShape& Move(const Vector& vector) override;

  bool ContainsPoint(const Point& point) const override;

  bool CrossesSegment(const Segment& segment) const override;

  IShape* Clone() const override;

  std::string ToString() override;

 private:
  Point start_point_;
  Vector direction_vector_;
};

class Line : public IShape {
 public:
  Line();
  Line(const Line& line);
  Line(const Point& first_point, const Point& second_point);
  Line(const int& a, const int& b, const int& c);

  IShape& Move(const Vector& vector) override;

  bool ContainsPoint(const Point& point) const override;

  bool CrossesSegment(const Segment& segment) const override;

  IShape* Clone() const override;

  std::string ToString() override;

  double DistanceFromPoint(const Point& point) const;

 private:
  int a_;
  int b_;
  int c_;
  Point first_point_;
  Point second_point_;
};

class Polygon : public IShape {
 public:
  explicit Polygon(const std::vector<Point>& points);

  IShape& Move(const Vector& vector) override;

  bool ContainsPoint(const Point& point) const override;

  bool CrossesSegment(const Segment& segment) const override;

  IShape* Clone() const override;

  std::string ToString() override;

 private:
  std::vector<Point> points_;
};

class Circle : public IShape {
 public:
  Circle(const Point& point, const int& radius);

  IShape& Move(const Vector& vector) override;

  bool ContainsPoint(const Point& point) const override;

  bool CrossesSegment(const Segment& segment) const override;

  IShape* Clone() const override;

  std::string ToString() override;

 private:
  Point center_point_;
  int radius_;
};

/*-----------------------------------Vector-----------------------------------*/

Vector::Vector() {
  x = 0;
  y = 0;
}

Vector::Vector(const CoordinateType& p_x, const CoordinateType& p_y) {
  x = p_x;
  y = p_y;
}

Vector::Vector(const Vector& vector) {
  x = vector.x;
  y = vector.y;
}

Vector::Vector(const Point& first_point, const Point& second_point) {
  x = first_point.point_vector.x - second_point.point_vector.x;
  y = first_point.point_vector.y - second_point.point_vector.y;
}

Vector Vector::operator+=(const Vector& vector) {
  this->x += vector.x;
  this->y += vector.y;
  return *this;
}

Vector Vector::operator-=(const Vector& vector) {
  this->x -= vector.x;
  this->y -= vector.y;
  return *this;
}

Vector Vector::operator*=(const MultiplierType& mul) {
  this->x *= mul;
  this->y *= mul;
  return *this;
}

Vector Vector::operator-() {
  this->x = -this->x;
  this->y = -this->y;
  return *this;
}

Vector& Vector::operator=(const Vector& vector) {
  if (this == &vector) {
    return *this;
  }
  x = vector.x;
  y = vector.y;
  return *this;
}

Vector::OperationResultType operator*(const Vector& v1, const Vector& v2) {
  return v1.x * v2.x + v1.y * v2.y;
}

Vector operator*(const Vector& v1, const Vector::MultiplierType& mul) {
  Vector vector(v1.x * mul, v1.y * mul);
  return vector;
}

Vector operator*(const Vector::MultiplierType& mul, const Vector& v1) {
  return v1 * mul;
}

Vector::OperationResultType operator^(const Vector& v1, const Vector& v2) {
  return v1.x * v2.y - v2.x * v1.y;
}

Vector operator+(const Vector& v1, const Vector& v2) {
  Vector vector(v1.x + v2.x, v1.y + v2.y);
  return vector;
}

Vector operator-(const Vector& v1, const Vector& v2) {
  Vector vector(v1.x - v2.x, v1.y - v2.y);
  return vector;
}

Vector::OperationResultType Vector::GetScalarSquare() const {
  return (*this) * (*this);
}

/*-----------------------------------Point------------------------------------*/

Point::Point(const Vector& vector) {
  point_vector.x = vector.x;
  point_vector.y = vector.y;
}

Point::Point(const Vector::CoordinateType& p_x,
             const Vector::CoordinateType& p_y) {
  point_vector.x = p_x;
  point_vector.y = p_y;
}

Vector operator-(const Point& first_point, const Point& second_point) {
  return {first_point.point_vector.x - second_point.point_vector.x,
          first_point.point_vector.y - second_point.point_vector.y};
}

IShape& Point::Move(const Vector& vector) {
  point_vector.x += vector.x;
  point_vector.y += vector.y;
  return *this;
}

bool Point::ContainsPoint(const Point& point) const {
  return (point_vector.x == point.point_vector.x &&
      point_vector.y == point.point_vector.y);
}

bool Point::CrossesSegment(const Segment& segment) const {
  return segment.ContainsPoint(*this);
}

IShape* Point::Clone() const {
  auto* point = new Point(point_vector);
  return point;
}

std::string Point::ToString() {
  return "Point(" + std::to_string(point_vector.x) + ", " +
      std::to_string(point_vector.y) + ")";
}

Point& Point::operator=(const Point& point) {
  if (this == &point) {
    return *this;
  }
  point_vector = point.point_vector;
  return *this;
}

bool operator==(const Point& first_point, const Point& second_point) {
  return (first_point.point_vector.x == second_point.point_vector.x) &&
      (first_point.point_vector.y == second_point.point_vector.y);
}

bool operator!=(const Point& first_point, const Point& second_point) {
  return (first_point.point_vector.x != second_point.point_vector.x) ||
      (first_point.point_vector.y != second_point.point_vector.y);
}

/*----------------------------------Segment-----------------------------------*/

Segment::Segment(const Point& first_point, const Point& second_point) {
  first_point_ = first_point;
  second_point_ = second_point;
}

Segment::Segment(const Segment& segment) {
  first_point_ = segment.first_point_;
  second_point_ = segment.second_point_;
}

IShape& Segment::Move(const Vector& vector) {
  first_point_.Move(vector);
  second_point_.Move(vector);
  return *this;
}

bool Segment::ContainsPoint(const Point& point) const {
  Vector ax_vector(first_point_, point);
  Vector bx_vector(second_point_, point);

  return !(bool)(ax_vector ^ bx_vector) &&
      point.point_vector.x <= std::max(first_point_.point_vector.x,
                                       second_point_.point_vector.x) &&
      point.point_vector.x >= std::min(first_point_.point_vector.x,
                                       second_point_.point_vector.x) &&
      point.point_vector.y <= std::max(first_point_.point_vector.y,
                                       second_point_.point_vector.y) &&
      point.point_vector.y >= std::min(first_point_.point_vector.y,
                                       second_point_.point_vector.y);
}

bool Segment::CrossesSegment(const Segment& segment) const {
  auto points = segment.GetPoints();

  Vector::OperationResultType first_mul =
      (GetDirectionVector() ^ Vector(GetPoints().first, points.second)) *
          (GetDirectionVector() ^ Vector(GetPoints().first, points.first));
  Vector::OperationResultType second_mul =
      (segment.GetDirectionVector() ^ Vector(points.first, second_point_)) *
          (segment.GetDirectionVector() ^ Vector(points.first, first_point_));

  return (first_mul < 0 && second_mul < 0) ||
      ContainsPoint(segment.GetPoints().first) ||
      ContainsPoint(segment.GetPoints().second) ||
      segment.ContainsPoint(GetPoints().first) ||
      segment.ContainsPoint(GetPoints().second);
}

IShape* Segment::Clone() const {
  auto* segment = new Segment(*this);
  return segment;
}
std::string Segment::ToString() {
  return "Segment(" + first_point_.ToString() + ", " +
      second_point_.ToString() + ")";
}

std::pair<Point, Point> Segment::GetPoints() const {
  return {first_point_, second_point_};
}

Vector Segment::GetDirectionVector() const {
  return second_point_ - first_point_;
}

/*------------------------------------Ray-------------------------------------*/

Ray::Ray(const Ray& ray) {
  start_point_ = ray.start_point_;
  direction_vector_ = ray.direction_vector_;
}

Ray::Ray(const Point& start_point, const Vector& direction_vector) {
  start_point_ = start_point;
  direction_vector_ = direction_vector;
}

Ray::Ray(const Point& first_point, const Point& second_point) {
  start_point_ = first_point;
  direction_vector_ = second_point - first_point;
}

IShape& Ray::Move(const Vector& vector) {
  start_point_.Move(vector);
  return *this;
}

bool Ray::ContainsPoint(const Point& point) const {
  return (direction_vector_ ^ (point - start_point_)) == 0 &&
      (direction_vector_ * (point - start_point_)) >= 0;
}

bool Ray::CrossesSegment(const Segment& segment) const {
  auto points = segment.GetPoints();
  if (ContainsPoint(points.first) || ContainsPoint(points.second)) {
    return true;
  }
  Vector first_vector = points.first - start_point_;
  Vector second_vector = points.second - start_point_;
  Line line(start_point_, points.first);
  Point c(start_point_.point_vector.x + direction_vector_.x,
          start_point_.point_vector.y + direction_vector_.y);

  return ((first_vector ^ direction_vector_) *
      (second_vector ^ direction_vector_) <
      0) &&
      (!line.CrossesSegment(Segment(points.second, c)));
}

IShape* Ray::Clone() const {
  auto* ray = new Ray(*this);
  return ray;
}

std::string Ray::ToString() {
  return "Ray(" + start_point_.ToString() + ", " + "Vector(" +
      std::to_string(direction_vector_.x) + ", " +
      std::to_string(direction_vector_.y) + "))";
}

/*------------------------------------Line------------------------------------*/

Line::Line() {
  a_ = 0;
  b_ = 0;
  c_ = 0;
  first_point_ = Point(Vector(0, 0));
  second_point_ = Point(Vector(0, 0));
}

Line::Line(const Line& line) {
  a_ = line.a_;
  b_ = line.b_;
  c_ = line.c_;
  first_point_ = line.first_point_;
  second_point_ = line.second_point_;
}

Line::Line(const Point& first_point, const Point& second_point) {
  a_ = second_point.point_vector.y - first_point.point_vector.y;
  b_ = first_point.point_vector.x - second_point.point_vector.x;
  c_ = first_point.point_vector.y *
      (second_point.point_vector.x - first_point.point_vector.x) -
      first_point.point_vector.x *
          (second_point.point_vector.y - first_point.point_vector.y);
  first_point_ = first_point;
  second_point_ = second_point;
}

Line::Line(const int& a, const int& b, const int& c) {
  a_ = a;
  b_ = b;
  c_ = c;
  first_point_ = Point(Vector(0, 0));
  second_point_ = Point(Vector(0, 0));
}

IShape& Line::Move(const Vector& vector) {
  c_ -= (a_ * vector.x + b_ * vector.y);
  return *this;
}

bool Line::ContainsPoint(const Point& point) const {
  return !(bool)(a_ * point.point_vector.x + b_ * point.point_vector.y + c_);
}

bool Line::CrossesSegment(const Segment& segment) const {
  return ((second_point_ - first_point_) ^
      (segment.GetPoints().first - first_point_)) *
      ((second_point_ - first_point_) ^
          (segment.GetPoints().second - first_point_)) <=
      0;
}

IShape* Line::Clone() const {
  auto* line = new Line(*this);
  return line;
}

std::string Line::ToString() {
  return "Line(" + std::to_string(a_) + ", " + std::to_string(b_) + ", " +
      std::to_string(c_) + ")";
}

double Line::DistanceFromPoint(const Point& point) const {
  return abs(a_ * point.point_vector.x + b_ * point.point_vector.y + c_) /
      sqrt(a_ * a_ + b_ * b_);
}

/*----------------------------------Polygon-----------------------------------*/

Polygon::Polygon(const std::vector<Point>& points) { points_ = points; }

IShape& Polygon::Move(const Vector& vector) {
  for (auto& point : points_) {
    point.Move(vector);
  }
  return *this;
}

bool Polygon::ContainsPoint(const Point& point) const {
  int counter;
  Ray ray;
  int n = (int)points_.size();

  do {
    ray = Ray(point, Vector((rand() & 666) + 1, (rand() & 666) + 1));
    counter = 0;
    for (const auto& point_of : points_) {
      if (point == point_of) {
        return true;
      }
      if (!ray.ContainsPoint(point_of)) {
        counter++;
      }
    }
  } while (counter != n);

  counter = 0;

  for (int i = 0; i < n; i++) {
    if (i == n - 1) {
      if (Segment(points_[i], points_[0]).ContainsPoint(point)) {
        return true;
      }
      if (ray.CrossesSegment(Segment(points_[i], points_[0]))) {
        counter++;
      }
      break;
    }
    if (Segment(points_[i], points_[i + 1]).ContainsPoint(point)) {
      return true;
    }
    if (ray.CrossesSegment(Segment(points_[i], points_[i + 1]))) {
      counter++;
    }
  }
  return (bool)(counter % 2);
}

bool Polygon::CrossesSegment(const Segment& segment) const {
  for (size_t i = 0; i < points_.size(); i++) {
    if (i == points_.size() - 1) {
      if (Segment(points_[i], points_[0]).CrossesSegment(segment)) {
        return true;
      }
      break;
    }
    if (Segment(points_[i], points_[i + 1]).CrossesSegment(segment)) {
      return true;
    }
  }
  return false;
}

IShape* Polygon::Clone() const {
  auto* polygon = new Polygon(*this);
  return polygon;
}

std::string Polygon::ToString() {
  std::string string = "Polygon(";
  for (size_t i = 0; i < points_.size(); i++) {
    string += points_[i].ToString();
    if (i == points_.size() - 1) {
      string += ")";
      break;
    }
    string += ", ";
  }
  return string;
}

/*-----------------------------------Circle-----------------------------------*/

Circle::Circle(const Point& point, const int& radius) {
  center_point_ = point;
  radius_ = radius;
}

IShape& Circle::Move(const Vector& vector) {
  center_point_.Move(vector);
  return *this;
}

bool Circle::ContainsPoint(const Point& point) const {
  return (point.point_vector.x - center_point_.point_vector.x) *
      (point.point_vector.x - center_point_.point_vector.x) +
      (point.point_vector.y - center_point_.point_vector.y) *
          (point.point_vector.y - center_point_.point_vector.y) <=
      radius_ * radius_;
}

bool Circle::CrossesSegment(const Segment& segment) const {
  Vector first_vector = segment.GetPoints().first - center_point_;
  Vector second_vector = segment.GetPoints().second - center_point_;
  Vector::OperationResultType first_scalar_square =
      std::min(first_vector.GetScalarSquare(), second_vector.GetScalarSquare());
  Vector::OperationResultType second_scalar_square =
      std::max(first_vector.GetScalarSquare(), second_vector.GetScalarSquare());

  if (first_scalar_square <= radius_ * radius_ &&
      second_scalar_square >= radius_ * radius_) {
    return true;
  }

  if (first_scalar_square >= radius_ * radius_ &&
      second_scalar_square >= radius_ * radius_) {
    Line line(segment.GetPoints().first, segment.GetPoints().second);
    return line.DistanceFromPoint(center_point_) <= radius_;
  }
  return false;
}

IShape* Circle::Clone() const {
  auto* circle = new Circle(*this);
  return circle;
}

std::string Circle::ToString() {
  return "Circle(" + center_point_.ToString() + ", " + std::to_string(radius_) +
      ")";
}

/*----------------------------------THE_END-----------------------------------*/
}  // namespace Geometry
