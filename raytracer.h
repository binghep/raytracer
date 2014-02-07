


struct Vec3  //vector3D
{
	float x, y, z;
	Vec3(float X = 0, float Y = 0, float Z = 0)
	{
		x = X;
		y = Y;
		z = Z;
	}
	float GetMagnitude()
	{
		return sqrt(x * x + y * y + z * z);
	}
	Vec3 operator*(float num) const//vec3*float=vec3
	{
		return Vec3(x * num, y * num, z * num);
	}
	Vec3 operator/(float num) const//vec3/float=vec3
	{
		return Vec3(x / num, y / num, z / num);
	}
	Vec3 operator +(const Vec3& v) const//vec3+vec3=vec3
	{
		return Vec3(x + v.x, y + v.y, z + v.z);
	}
	Vec3 operator -(const Vec3& v) const//vec3-vec3=vec3
	{
		return Vec3(x - v.x, y - v.y, z - v.z);
	}
	//operator overloading
	//cross product using %
	Vec3 operator %(const Vec3& v) const
	{
		return Vec3(y * v.z - z * v.y,
			z * v.x - x * v.z,
			x * v.y - y * v.x);
	}
	//dot product using &
	float operator &(const Vec3& v) const
	{
		return (x * v.x + y * v.y + z * v.z);
	}
};


struct Point{
	float x, y, z;
	Point(float X = 0, float Y = 0, float Z = 0)
	{
		x = X;
		y = Y;
		z = Z;
	}
	Vec3 operator -(const Point& v) const//point-point=vec3
	{
		return Vec3(x - v.x, y - v.y, z - v.z);
	}
	Point operator +(const Point& v) const//point+point=point
	{
		return Point(x + v.x, y + v.y, z + v.z);
	}
	Point operator +(const Vec3& v) const//point+vec3=point
	{
		return Point(x + v.x, y + v.y, z + v.z);
	}
	Point operator -(const Vec3& v) const//point-vec3=point
	{
		return Point(x - v.x, y - v.y, z - v.z);
	}
	Point operator /(const float& num) const//point/float=point
	{
		return Point(x / num, y / num, z / num);
	}
	Point operator *(const float& num) const//point*float=point
	{
		return Point(x * num, y * num, z * num);
	}
	/*double equal operator*/
	bool operator== (const Point &c1) {
		return (x == c1.x && y == c1.y && z == c1.z);
	}
};

struct ColorType
{
	float r, g, b;
	ColorType(float X = 0, float Y = 0, float Z = 0)
	{
		r = X;
		g = Y;
		b = Z;
	}
	bool operator== (const ColorType &c1) {
		return (r == c1.r && g == c1.g && b == c1.b);
	}
	bool operator!= (const ColorType &c1) {
		return (r != c1.r || g != c1.g || b != c1.b);
	}
	/*ColorType & operator +=(const ColorType& rhs)
	{
	// actual addition of rhs to *this
	this->r = this->r + rhs.r;
	this->g = this->g + rhs.g;
	this->b = this->b + rhs.b;
	return *this;
	}*/
	ColorType operator +(const ColorType& c1) const//get new ColorType
	{
		return ColorType(r + c1.r, g + c1.g, b + c1.b);
	}
	ColorType operator *(const float& num) const//get new ColorType
	{
		return ColorType(r*num, g*num, b*num);
	}

};

struct RayType {
	float x, y, z;
	float dx, dy, dz;

};


struct SphereType {

	float x, y, z, r;
	Point center;
	//float x1,y1,z1;//intersection with a ray

};

struct LightType{
	float x, y, z, r, g, b;//r g b are color intensity
	//usually 0.5 each
	int w;//0 for directional light, 1 for point light
	Point lightlocation;
	Vec3 directionallight;
};
struct MaterialColor{
	float dr, dg, db, sr, sg, sb, ka, kd, ks, n, alpha, eta;
};
struct Node{
	int n;
	float p1, p2;
};
struct s{
	RayType parent, ch, ch2;
	float afterin, beforeout, afterout, Frout,transmitconstant;
	Point in, out;
	Vec3 N, N1, T,d,d2;
	s(RayType* X){
		parent = *X;
	}
	void print();
	bool child(int whichsphere);
};
float shadow(RayType *ray, float distance,int whichsphere);
ColorType Trace_Ray(RayType ray, Vec3 ray_dir, Point from);
#define MAXLIGHT 4
#define MAXSPHERE 20
ColorType reflective_trace_and_shade(RayType *ray);
void pointcolor(Point *p);
void add_two_colors_due_to_this_ray_in_sphere(RayType *trans,int whichsphere);
