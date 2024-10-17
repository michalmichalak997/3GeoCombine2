#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <vector>


using namespace std;


const int n = 3;
const int g = 3;
const double coeff = 180 / M_PI;
vector <double> z_axis={0,0,1};

int n_comb(int s)
{
	int komb = s;
	for (int i = s; i > s - 2; i--)
	{
		komb = komb * (i - 1);
	}

	return komb / 6;
}

class Point
{
    private:
	double x, y, z;


	public:
    Point() : x(-9999), y(-9999), z(-9999) {} // Default constructor

	Point(double x_1, double x_2, double x_3)
		:x(x_1), y(x_2), z(x_3) {}

    double get_x() const //a
	{
return this->x;
	}

	double get_y() const //a
	{
return this->y;
	}

	double get_z() const //a
	{
return this->z;
	}

};

class plane //a class that stores the crucial figures in terms of computing the orientation
{


private:
    Point first_point;
    Point second_point;
    Point third_point;
    vector <double> first_vec;
    vector <double> second_vec;
    vector <double> third_vec;
    vector <double> normal_vec;			//normal vector of a triangle
    vector <double> directional;
    vector <double> dip_vec;
    double doc;
    bool lin_dependence;
    string dip_degrees;
    string azimuth_degrees;

public:
	const double ex = 2;			//we introduce the restriction of collinearity
	//double first_vec[n];            //the first edge of a triangle
	//double second_vec[n];			//the second edge of a triangle
	//double third_vec[n];			//the third edge of a triangle

	//double directional[n];			//the projection of the normal vector onto the horizontal plane
	//double z_axis[n] = { 0,0,1 };   //the definition of the z-axis
	//double doc;						//a variable that contains the collinearity coefficient
	//double area;					//a variable that stores the area of a triangle
	//bool lin_dependence;		    //a bool variable to check to answer whether points are (too) collinear
	//string dip_degrees;             //a text variable to store the dip angle
	//string azimuth_degrees;         //a text variable to store the dip direction

	double dip_azimuth(vector <double> normal, int n = 2) //a function that computes the dip azimuth
	{

		double angle = atan2(normal[1], normal[0]);
		angle = angle * coeff;
		if (angle < 0)
		{
			return angle + 360;
		}
		else
		{
			return angle;
		}
	}

	double dip_angle(vector <double> z_axis, vector <double> normal_v) //function that computes the dip angle
	{
		double numerator=dot_product(z_axis, normal_v);
		if (numerator>=1.0){
            return 0;
		}
		else {
            return acos(numerator)*coeff;
		}
	}

	static double dot_product(vector <double> vector_1, vector <double> vector_2, int n = 3) //function that computes the dot product of vectors
	{
		double product = 0;
		for (int i = 0; i < n; i++)
		{
			product += vector_1[i] * vector_2[i];
		}
		return product;
	}

	bool dependence(vector <double> v1, vector <double> v2, vector <double> v3) //function that checks whether the points are collinear
	{
		double len_v1 = length(v1);
		double len_v2 = length(v2);
		double len_v3 = length(v3);
		double lengths[n] = { len_v1, len_v2, len_v3 };

		sort(lengths, lengths + n);
		this->doc = lengths[2] / (lengths[0] + lengths[1]);
		int k = 0;
		for (int i = 0; i < n; i++)
		{
			if (lengths[i] == 0)
			{
			    k+=1;
				//throw runtime_error("Points coincidence");
			}
		}
		if (k != 0)
		{
			return true;
		}
		else
		{
			if (doc > ex)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	static double length(vector <double> line_vector, int n = 3) //function that computes the length of a vector
	{
		double vector_length = sqrt(pow(line_vector[0], 2) + pow(line_vector[1], 2) + pow(line_vector[2], 2));
		return vector_length;
	}

    vector <string> center() //function that computes the geometrical centre of a triangle, usuwane
	{
		double x = (first_point.get_x() + second_point.get_x() + third_point.get_x()) / (3.0);
		double y = (first_point.get_y() + second_point.get_y() + third_point.get_y()) / (3.0);
		double z = (first_point.get_z() + second_point.get_z() + third_point.get_z()) / (3.0);
		vector<string> napis{ to_string(x), to_string(y), to_string(z) };
		return napis;
	}

	plane( Point point_1, Point point_2, Point point_3) //the class constructor
	{
	    this->first_point=Point(point_1.get_x(), point_1.get_y(), point_1.get_z());
		this->second_point=Point(point_2.get_x(), point_2.get_y(), point_2.get_z());
		this->third_point=Point(point_3.get_x(), point_3.get_y(), point_3.get_z());

		vector<double>  first_try = { point_2.get_x() - point_1.get_x(), point_2.get_y() - point_1.get_y(), point_2.get_z() - point_1.get_z() };
		vector<double>  second_try = { point_3.get_x() - point_1.get_x(), point_3.get_y() - point_1.get_y(), point_3.get_z() - point_1.get_z() };
		vector<double>  third_try = { point_3.get_x() - point_2.get_x(), point_3.get_y() - point_2.get_y(), point_3.get_z() - point_2.get_z() };


		bool test = dependence(first_try, second_try, third_try);
		if (test == true)
		{
			lin_dependence = true;
			this->first_vec = { -99,-99,-99 };
			this->second_vec = { -99,-99,-99 };
			this->third_vec = { -99,-99,-99 };
			this->normal_vec = { -99,-99,-99 };
			this->dip_vec = { -99, -99, -99 };
			measure();
		}
		else
		{
			lin_dependence = false;
			for (int i = 0; i < n; i++)
			{
				this->first_vec.push_back(first_try[i]);
				this->second_vec.push_back(second_try[i]);
				this->third_vec.push_back(third_try[i]);
			}

			normal_vec.push_back(first_vec[1] * second_vec[2] - second_vec[1] * first_vec[2]);
			normal_vec.push_back(first_vec[2] * second_vec[0] - second_vec[2] * first_vec[0]);
			normal_vec.push_back(first_vec[0] * second_vec[1] - second_vec[0] * first_vec[1]);
			if (normal_vec[2] < 0) {
				normal_vec[0] *= -1;
				normal_vec[1] *= -1;
				normal_vec[2] *= -1;
			}

			double normal_length=length(normal_vec);

			normal_vec[0] /= normal_length;
			normal_vec[1] /= normal_length;
			normal_vec[2] /= normal_length;

			if (normal_vec[2] == 1) {
				dip_vec.push_back(-9999);
				dip_vec.push_back(-9999);
				dip_vec.push_back(-9999);
			}
			else {
				dip_vec.push_back(cos(dip_angle(z_axis, normal_vec) / coeff) * cos(dip_azimuth(normal_vec) / coeff));
				dip_vec.push_back(cos(dip_angle(z_axis, normal_vec) / coeff) * sin(dip_azimuth(normal_vec) / coeff));
				dip_vec.push_back(-sin(dip_angle(z_axis, normal_vec) / coeff));
			}

			measure(); //setting dip_degrees and azimuth_degrees


		}


	}


	void measure()//function that supplies orientation results also for singularities
	{
		if (lin_dependence) // result for collinear points
		{
			this->azimuth_degrees = ("LT");
			this->dip_degrees = ("LT");

		}
		else if (normal_vec[0] == 0 && normal_vec[1] == 0 && normal_vec[2] != 0) //result for a horizontal triangle
		{
			this->dip_degrees = "0";
			this->azimuth_degrees = ("x");

		}
		else if (normal_vec[2] == 0) //result for a vertical triangle
		{
			this->dip_degrees = "90";
			this->azimuth_degrees = to_string(dip_azimuth(normal_vec));
		}
		else //a normal case (no singularities)
		{
			double dipping_angle = dip_angle(z_axis, normal_vec);
			this->dip_degrees = to_string(dipping_angle);
			this->azimuth_degrees = to_string(dip_azimuth(normal_vec));
		}
	}


	string get_dip_angle() {

		return this->dip_degrees;
	}

	string get_azimuth() {

		return this->azimuth_degrees;
	}

	vector<double> get_normal() //function that -computes- returns the normal vector
	{
		vector<double> normal_vector = { this->normal_vec[0] ,this->normal_vec[1], this->normal_vec[2] };

		return normal_vector;
	}

	double get_doc()
	{
		return this->doc;
	}

	Point get_first_point()
	{
		return this->first_point;
	}

	Point get_second_point()
	{
		return this->second_point;
	}

	Point get_third_point()
	{
		return this->third_point;
	}

	vector<double> get_first_vec()
	{
		return this->first_vec;
	}

	vector<double> get_second_vec()
	{
		return this->second_vec;
	}
	vector<double> get_third_vec()
	{
		return this->third_vec;
	}

	vector<double> get_directional()
	{
		return this->directional;
	}

	vector<double> get_dip() {

		vector<double> dip_vector = { dip_vec[0],dip_vec[1] ,dip_vec[2] };
		return(dip_vector);
	}

	bool get_lin_dependence() {

		return this->lin_dependence;
	}






};

int main()
{
		try
		{
			string path_i, path_o;

			cout << "Type in the path of your input data:" << endl; //the user is required to type in the input path
			cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\JurassicBottomInput.txt" << endl << endl;
			cin >> path_i;

			ifstream download(path_i);
			if (!download) cout << "Error in opening file" << endl; //the case when the file cannot be uploaded

			cout << "Type in the path of the output:" << endl; //the user is required to type in the output path
			cout << "Example: C:\\dev\\CGAL-4.8\\examples\\Triangulation_2\\JurassicBottomOutput.txt" << endl << endl;
			cin >> path_o;

			string tempor;//a temporary variable storing figures while uploading
			stringstream ss;
            vector <std::pair< Point, int> > pts; //pits replaced by pts
            int number_of_point=0;
			while (getline(download, tempor))//loading points line-by-line
			{
				istringstream convert(tempor);
				double a, b, c;
				if (!(convert >> a >> b >> c)) { cerr << "Error downloading data - check!"; }
				pts.push_back(  std::make_pair(Point(a, b, c), number_of_point+1)  );
				number_of_point+=1;
			}

			int n_boreholes= pts.size();
			int l_komb = n_comb(n_boreholes);

			if (n_boreholes < 3)  throw runtime_error("You need at least three boreholes");

			cout << "You have successfully uploaded " << n_boreholes << " boreholes" <<  endl;
			cout << "This will give: " << l_komb << " observations." << endl;

			ofstream all(path_o);

			all << "X1" << ';' << "Y1" << ';' << "Z1" << ';' << "X2" << ';' <<
				"Y2" << ';' << "Z2" << ';' << "X3" << ';' << "Y3" << ';' << "Z3"
				<< ';' << "X_C" << ';' << "Y_C" << ';' << "Z_C" << ';' << "X_N" << ';' << "Y_N"
				<< ';' << "Z_N" << ';' << "Dip_ang" << ';' << "Dip_dir" <<
				';' << "DOC" << ';' << "ID1" << ';' << "ID2" <<  ';' << "ID3" <<  '\n';

			if (n_boreholes == 3) {


				Point point_1(pts.at(0).first.get_x(), pts.at(0).first.get_y(), pts.at(0).first.get_z());
				Point point_2(pts.at(1).first.get_x(), pts.at(1).first.get_y(), pts.at(1).first.get_z());
				Point point_3(pts.at(2).first.get_x(), pts.at(2).first.get_y(), pts.at(2).first.get_z());

				plane current_plane = plane(point_1, point_2, point_3);					//constructing a plane that is processed at the moment
				//string result = current_plane.measure();								//extracting the dip angle and the dip direction
				vector<string> centroid = current_plane.center();


				all << to_string(point_1.get_x()) << ";" << to_string(point_1.get_y() ) << ";" << to_string(point_1.get_z()) << ";" << //saving orientation elements with respect to the column names
					to_string(point_2.get_x() ) << ";" << to_string(point_2.get_y() ) << ";" << to_string(point_2.get_z() ) << ";" <<
					to_string(point_3.get_x() ) << ";" << to_string(point_3.get_y() ) << ";" << to_string(point_3.get_z() ) << ";" <<
					centroid[0] << ";" << centroid[1] << ";" << centroid[2] << ";" << current_plane.get_normal()[0] << ";" << current_plane.get_normal()[1] << ";" << current_plane.get_normal()[2] <<
					 ";" << current_plane.get_dip_angle()  << ";" << current_plane.get_azimuth() << ";" << current_plane.get_doc() << endl;

			}
			else

			{
					int S[g];
					for (int i = 0; i < g; i++) S[i] = i;

					int p = g - 1;

					while (p >= 0)
					{
						Point point_1(pts.at(S[0]).first.get_x(), pts.at(S[0]).first.get_y(), pts.at(S[0]).first.get_z());
                        Point point_2(pts.at(S[1]).first.get_x(), pts.at(S[1]).first.get_y(), pts.at(S[1]).first.get_z());
                        Point point_3(pts.at(S[2]).first.get_x(), pts.at(S[2]).first.get_y(), pts.at(S[2]).first.get_z());

						plane current_plane = plane(point_1, point_2, point_3);					//constructing a plane that is processed at the moment
						//string result = current_plane.measure();								//extracting the dip angle and the dip direction
						vector<string> centroid = current_plane.center();



                    all << to_string(point_1.get_x()) << ";" << to_string(point_1.get_y() ) << ";" << to_string(point_1.get_z()) << ";" << //saving orientation elements with respect to the column names
					to_string(point_2.get_x() ) << ";" << to_string(point_2.get_y() ) << ";" << to_string(point_2.get_z() ) << ";" <<
					to_string(point_3.get_x() ) << ";" << to_string(point_3.get_y() ) << ";" << to_string(point_3.get_z()) << ";" <<
					centroid[0] << ";" << centroid[1] << ";" << centroid[2] << ";" << current_plane.get_normal()[0] << ";" << current_plane.get_normal()[1] << ";" << current_plane.get_normal()[2] <<
					 ";" << current_plane.get_dip_angle()  << ";" << current_plane.get_azimuth() << ";" << current_plane.get_doc() <<  ";" <<
					 to_string( pts.at(S[0]).second ) << ";" << to_string( pts.at(S[1]).second) << ";" << to_string(pts.at(S[2]).second )  << endl;

						if (S[g - 1] == n_boreholes - 1) { p--; }
						else { p = g - 1; }
						if (p >= 0)
						{
							for (int i = g - 1; i >= p; i--)
							{
								S[i] = S[p] + i - p + 1;
							}
						}

					}

			}

		}
		catch (runtime_error e)
		{

			cout << "Runtime error: " << e.what();
		}

	system("PAUSE");
	return 0;
}
