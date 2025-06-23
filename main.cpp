
/// Created by Karol Kubek
#include <bits/stdc++.h>
#define PI (2 * asinl(1))
#define MP make_pair
#define SQR(x) ((x)*(x))
#define fi first
#define se second

#define EPS 1e-7
#define FOR(v,p,k) for(int v=p;v<=k;++v)
#define FORD(v,p,k) for(int v=p;v>=k;--v)
#define REP(i,n) for(int i=0;i<(n);++i)
#define VAR(v,i) __typeof(i) v=(i)
#define PB push_back
using std::cout, std::cin, std::string, std::pair, std::vector, std::swap;
typedef long double ld;
typedef string st;
typedef pair<int,int> PII;
typedef vector<int> VI;
typedef struct {
  ld x;
  ld y;
} Point;
typedef struct {
 Point A;
 Point B;
} seg;

inline Point operator- (const Point &B, const Point &A) {
 Point res;
 res.x = A.x - B.x;
 res.y = A.y - B.y;
 return res;
}
inline bool operator ==(const Point &A, const Point &B) {
 if (A.x != B.x || A.y != B.y) return 0;
 return 1;
}

inline Point operator* (const ld lambda,const Point A)
{
    Point res;
    res.x = lambda*A.x;
    res.y = lambda*A.y;
}



inline int iszero(ld x) {
 return (x < EPS) && (x > -EPS);
}
const int number = 256;
char *tr(char *str) /// Funkcja pozwalaj¹ca mi u¿ywaæ polskich znaków w konsoli
  {                 /// Zapo¿yczona z forum Miros³awa Zelenta
   static char buff[number];
   char cp[]="\176\245\206\251\210\344\242\230\276\253\244\217\250\235\343\340\227\275\215°¹æê³ñóœ¿Ÿ¥ÆÊ£ÑÓŒ¯";
   if(strlen(str)>=sizeof(buff)) return str;
   char *bf=buff;
   while(*str)
     {
      char *pos=strchr(cp+19,*str);
      *(bf++)=pos?*(pos-19):*str;
      ++str;
     }
   *bf=0;
   return buff;
  }

struct Triple{
   ld p1; /// Zawsze bêdê uto¿samia³ z wierzcho³kiem A lub bokiem a lub katem alfa a tak¿e dowoln¹ kresk¹ wychodzac¹ z A (h_a, m_a, etc)
   ld p2; /// Zawsze bêdê uto¿samia³ z wierzcho³kiem B lub bokiem b lub katem beta a tak¿e dowoln¹ kresk¹ wychodzac¹ z B
   ld p3; /// Zawsze bêdê uto¿samia³ z wierzcho³kiem C lub bokiem c lub katem gamma a tak¿e dowoln¹ kresk¹ wychodzac¹ z C
};
inline Point Sum (Point B, Point A) {
 Point res;
 res.x = A.x + B.x;
 res.y = A.y + B.y;
 return res;
}
inline ld LengthOfSegment(Point A,Point B)
{
    return sqrt((A.x-B.x)*(A.x-B.x)+(A.y-B.y)*(A.y-B.y));
}
/*double area(Point A,Point B, Point C)
{
     double result;
     result = (B.x - A.x)*(C.y-A.y);
     result-= (C.x - A.x)*(B.y-A.y);
     result = fabs(result);
     result/=2;
     return result;
}*/

struct Triangle{
   Point A;
   Point B;
   Point C;
   ld Segment_a = LengthOfSegment(B,C);
   ld Segment_b = LengthOfSegment(A,C);
   ld Segment_c = LengthOfSegment(A,B);
};
struct Circle{
   Point O;
   ld R;

};
inline ld Area(Triangle T) /// Funkcja zwracajaca pole trojkata T
{
    ld result;
    result = (T.B.x - T.A.x)*(T.C.y - T.A.y);
    result -=  (T.C.x - T.A.x)*(T.B.y-T.A.y);
    result = fabs(result); /// gwarantuje, ze pole wyjdzie dodatnie, nawet gdy
    /// policzony wyznacznik okazal sie ujemny;


    return 0.5*result;
}
/*bool isPointInsideTriangle(Triangle T,Point P){
    /// Jesli pola trojkatow ABP,ACP i BCP sa dodatnie
    /// Czyli jesli P nie lezy na zadnym z bokow
      if(Area(T1)>0)&&(Area(T2)>0)&&(Area(T3>0))
      {     /// To jesli suma pol trojkatow ABP,ACP i BCP jest rowna polu trojkata ABC
          if(Area(T1)+Area(T2)+Area(T3)==Area(T))
          return true; /// Zwroc prawde, gdyz oznacza to ze P lezy wewnatrz trojkata ABC
      }
      return false;
}*/
void isPointInsideTriangle(Triangle T,Triangle T1,Triangle T2,Triangle T3)
{
    if((Area(T1)==0)||(Area(T2)==0)||(Area(T3)==0))
        cout<<"Punkt P lezy na jednym z bokow trojkata ABC"<<'\n';
    else if((Area(T1)>0)&&(Area(T2)>0)&&(Area(T3)>0))
      {     /// To jesli suma pol trojkatow ABP,ACP i BCP jest rowna polu trojkata ABC
          if(Area(T1)+Area(T2)+Area(T3)==Area(T))
            cout<<"Punkt P lezy wewnatrz trojkata ABC"<<'\n';
       }
}
struct Line{   /// Struktura prostej - 2 punkty wyznaczaj¹ j¹ jednoznacznie
    Point P1;
    Point P2;
};

inline ld Perimeter(Triangle T) /// Funkcja obliczajaca L - obwod trojkata
{
    ld result;
    result = LengthOfSegment(T.A,T.B)+LengthOfSegment(T.A,T.C)+LengthOfSegment(T.C,T.B);
    return result;
}
inline ld Halfperimeter(Triangle T) /// Funkcja obliczajaca p - polowe obwodu trojkata p = (a+b+c)/2
{
    return 0.5*Perimeter(T);
}
inline ld RadiusOfCircumCircle(Triangle T) /// Funkcja obliczajaca R - promien okregu opisanego na trojkacie
{
    ld result;
    ld a = LengthOfSegment(T.C,T.B);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);

    result = a*b*c/4;
    result/=Area(T);
     return result;
}
inline ld RadiusOfInCircle(Triangle T) /// Funkcja obliczajaca r - promien okregu wpisanego w trojkat
{
     return Area(T)/Halfperimeter(T);
}
struct line{Point P1,P2;};
inline ld DistanceFromLine(line k,Point P)
{
    ld result;
    Triangle T;
    T.A = k.P1;
    T.B = k.P2;
    T.C = P;
    result = Area(T);
    result/= LengthOfSegment(k.P1,k.P2);
    return result;
}

inline Point Centroid(Triangle T) /// Funkcja zwraca wspó³rzêdne œrodka ciê¿koœci trójk¹ta
{
    Point G;
    G.x = (T.A.x+T.B.x+T.C.x)/3;
    G.y = (T.A.y+T.B.y+T.C.y)/3;

    return G;

}

bool AreCollineare(Triangle T) /// Funkcja rozstrzyga, czy 3 punkty wyznaczaj¹ trójk¹t, czy te¿ s¹ wspó³liniowe - wówczas trójk¹t jest zdegenerowany
{
    if(Area(T)==0) return true;
    return false;
}
st TypeOfTriangle(Triangle T){   /// Ze wzglêdu na boki
    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    st which_one;
   if((a==b)&&(b==c)) which_one = tr("równoboczny");
   else if((a==b)||(a==c)||(b==c)) which_one = tr("równoramienny");
   else which_one = tr("ró¿noboczny");
   return which_one;
}


/*inline ld Median(Triangle T,Point X){ /// Funkcja zwraca d³ugoœæ œrodkowej   NIE DZIA£¥

    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    if(X==T.C)
    return 0.5*sqrt(2*(a*a+b*b)-c*c);
    if(X==T.A)
    return 0.5*sqrt(2*(c*c+b*b)-a*a);
    if(X==T.B)
    return 0.5*sqrt(2*(a*a+c*c)-b*b);
}*/
inline Point MidP(Point A,Point B)
{
    Point M;
    M.x = 0.5*(A.x+B.x);
    M.y = 0.5*(A.y+B.y);
    return M;
}
/*inline Median(Triangle T,Point X){

  return LengthOfSegment(T.A,MidP(T.A,T.B));
}*/

Triple angels(Triangle T)
{
    Triple Katy;

    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
   ld Alpha = acos((b*b+c*c-a*a)/(2*b*c));
   ld Beta = acos((a*a+c*c-b*b)/(2*a*c));
   ld Gamma = acos((b*b+a*a-c*c)/(2*a*b));
   Katy.p1 = Alpha;
   Katy.p2 = Beta;
   Katy.p3 = Gamma;
   return Katy;

}
st TypeOfTriangle_Angles(Triangle T) /// Ze wzglêdu na k¹ty
{
    st which_one;

      ld alf = angels(T).p1;
      ld bet = angels(T).p2;
      ld gamm = angels(T).p3;
  /*ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    */
    if(alf == 0.5*PI || bet == 0.5*PI || gamm == 0.5*PI) which_one = tr("prostok¹tny");
    else if(alf < 0.5*PI && bet < 0.5*PI && gamm < 0.5*PI) which_one = tr("ostrok¹tny");
    else if(alf > 0.5*PI || bet > 0.5*PI || gamm > 0.5*PI) which_one = tr("rozwartok¹tny");
    return which_one;

}
Triple LengthOfAngelBisector(Triangle T)
{
    Triple Triple1;
    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    ld p = Halfperimeter(T);
    Triple1.p1 = (2/(b+c))*sqrt(p*(p-a)*b*c);
    Triple1.p2 = (2/(a+c))*sqrt(p*(p-b)*a*c);
    Triple1.p3 = (2/(a+b))*sqrt(p*(p-c)*a*b);
    return Triple1;
}
Point wektor(Point P1,Point P2)/// Zwraca wspó³rzêdne wektora P1P2 wektor (UWAGA! kolejnosc pkt istotna!)
{

    Point res;
    res.x = P2.x - P1.x;
    res.y = P2.y - P1.y;
    return res;


}

ld det(pair<ld,ld> a,pair<ld,ld> b) ///Prymitywna funkcja sadge oblicza wyznacznik 2x2
{
    return a.fi*b.se - a.se*b.fi;
}

Point Solve(pair<ld,ld> px,pair<ld,ld> py,pair<ld,ld> pconst) ///Zakladam, ze uklad jest oznaczony - no ale tak jest, bo trojkat istnieje xD
 {
     Point P;
     ld W = det(px,py);
     ld Wx = det(pconst,px);
     ld Wy = det(px,pconst);
     P.x = Wx/W;
     P.y = Wy/W;
     return P;
 }


Point FindCircumCenter(Triangle T) ///Dziala niepoprawnie!
{                                  /// Moja autorska - temu nie dzia³a xD
    Point P3,P4;
    Point O;
    ld x0 = MidP(T.A,T.B).x;
    ld y0 = MidP(T.A,T.B).y;
    P3 = wektor(T.A,T.B);
    ld a = P3.x;
    ld b = P3.y;

    ld x1 = MidP(T.A,T.C).x;
    ld y1 = MidP(T.A,T.C).y;

    P4 = wektor(T.A,T.C);
    ld c = P4.x;
    ld d = P4.y;

    pair<ld,ld> px = {a,c};
    pair<ld,ld> py = {b,d};
    ld w1 = a*x0+b*y0;
    ld w2 = c*x1+d*y1;

    pair<ld,ld> pconst = {w1,w2};


    O = Solve(px,py,pconst);
    return O;
}

Point CircleThroughThreePoints(Triangle T) /// Ta jest legit
{                                          /// Zapo¿yczona z geeksforgeeks.com
    Point O;
    ld x12 = T.A.x - T.B.x;
    ld x13 = T.A.x - T.C.x;

    ld y12 = T.A.y - T.B.y;
    ld y13 = T.A.y - T.C.y;

    ld y31 = T.C.y - T.A.y;
    ld y21 = T.B.y - T.A.y;

    ld x31 = T.C.x - T.A.x;
    ld x21 = T.B.x - T.A.x;

    ld sx13 = pow(T.A.x, 2) - pow(T.C.x, 2);

    // y1^2 - y3^2
    ld sy13 = pow(T.A.y, 2) - pow(T.C.y, 2);

    ld sx21 = pow(T.B.x, 2) - pow(T.A.x, 2);
    ld sy21 = pow(T.B.y, 2) - pow(T.A.y, 2);

    ld f = ((sx13) * (x12)
             + (sy13) * (x12)
             + (sx21) * (x13)
             + (sy21) * (x13))
            / (2 * ((y31) * (x12) - (y21) * (x13)));
    ld g = ((sx13) * (y12)
             + (sy13) * (y12)
             + (sx21) * (y13)
             + (sy21) * (y13))
            / (2 * ((x31) * (y12) - (x21) * (y13)));

    ld c = -pow(T.A.x, 2) - pow(T.A.y, 2) - 2 * g * T.A.x - 2 * f * T.A.y;

    // eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
    // where centre is (h = -g, k = -f) and radius r
    // as r^2 = h^2 + k^2 - c
    ld h = -g;
    ld k = -f;
    O.x = h;
    O.y = k;
   return O;
}
Point Incenter(Triangle T)
{
    Point res;
    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    ld p = a+b+c;
    res.x = (a*T.A.x+b*T.B.x+c*T.C.x)/p;
    res.y = (a*T.A.y+b*T.B.y+c*T.C.y)/p;
    return res;
}
Point Orthocenter(Triangle T)
{
    Point H;
    Point O = CircleThroughThreePoints(T);
    H.x = T.A.x+T.B.x+T.C.x-2*O.x;
    H.y = T.A.y+T.B.y+T.C.y-2*O.y;

    return H;

}
inline bool isPointInteger(Point P)
{
    if(iszero(P.x - floor(P.x))==true && iszero(P.y - floor(P.y))==true)
        return true;
}
bool areVerticesIntegerPoints(Triangle T)
{
    bool p = isPointInteger(T.A);
    bool q = isPointInteger(T.B);
    bool r = isPointInteger(T.C);
    if((p && q && r)==true)return true;
    return false;
}
/// Funckcja oblicza NWD dwoch liczb
int gcd(int a, int b)
{
    if((a==1)||(b==1))return 1;
    if(a%b==0)return b;
    if(b%a==0)return a;
    else{
    while(b)
    {
        swap(a%=b,b);
    }
    }
    return a;
}
int getBoundaryCount(Point P,Point Q)
{
    int p1 = (int)P.x;
    int p2 = (int)P.y;

    int q1 = (int)Q.x;
    int q2 = (int)Q.y;
    if (p1==q1)
        return abs(p2 - q2) - 1;
    if (p2 == q2)
        return abs(p1 - q1) - 1;

    return gcd(abs(p1-q1), abs(p2-q2)) - 1;
}
int getInternalCount(Triangle T)
{
    // 3 extra integer points for the vertices
    int BoundaryPoints = getBoundaryCount(T.A, T.B) +
                         getBoundaryCount(T.B, T.C) +
                         getBoundaryCount(T.A, T.C) + 3;

    // Calculate 2*A for the triangle
    int doubleArea = 2*Area(T);

    // Use Pick's theorem to calculate the no. of Interior points
    return (doubleArea - BoundaryPoints + 2)/2;

}
inline int integersPointsOnSides(Triangle T)
{
    int I = getInternalCount(T);
    int A = (int)Area(T);
    return (2*A+2-2*I);
}
Point Lemoine_Point(Triangle T)
{
    Point res;
    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);


    res.x = (a*a*T.A.x+b*b*T.B.x+c*c*T.C.x)/(a*a+b*b+c*c);
    res.y = (a*a*T.A.y+b*b*T.B.y+c*c*T.C.y)/(a*a+b*b+c*c);

    return res;
}
Point Gergonne_Point(Triangle T)
{
    Point res;
    ld s = Halfperimeter(T);
    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    res.x =((s-a)*(s-b)*T.A.x+(s-a)*(s-c)*T.B.x+(s-b)*(s-c)*T.C.x)/((s-a)*(s-b)+(s-a)*(s-c)+(s-b)*(s-c));
    res.y =((s-a)*(s-b)*T.A.y+(s-a)*(s-c)*T.B.y+(s-b)*(s-c)*T.C.y)/((s-a)*(s-b)+(s-a)*(s-c)+(s-b)*(s-c));

    return res;
}
Point Nagel_Point(Triangle T)
{
    Point res;
    ld s = Halfperimeter(T);
    ld a = LengthOfSegment(T.B,T.C);
    ld b = LengthOfSegment(T.A,T.C);
    ld c = LengthOfSegment(T.A,T.B);
    res.x = ((s-a)*T.A.x+(s-b)*T.B.x+(s-c)*T.C.x)/s;
    res.y = ((s-a)*T.A.y+(s-b)*T.B.y+(s-c)*T.C.y)/s;
    return res;
}
int main()
{

    std::ios_base::sync_with_stdio(0);

    Triangle T; /// Deklaruje istnienie mojej struktury - trójk¹ta
    cout<<"Created by Karol Kubek (C) 2022"<<'\n';
    cout<<tr("Niniejszy program dedykujê mojemu przyjacielowi Kacperkowi J. :)")<<'\n';
    cout<<'\n';
    cout<<tr("Podaj wspó³rzêdne wierzcho³ków trójk¹ta, oddzielone spacjami:")<<'\n';
    cout<<tr("Wspó³rzêdne pierwszego wierzcho³ka: ");
    cin >> T.A.x >> T.A.y;
    cout<<tr("Wspó³rzêdne drugiego wierzcho³ka: ");
    cin >> T.B.x >> T.B.y;
    cout<<tr("Wspó³rzêdne trzeciego wierzcho³ka: ");
    cin >> T.C.x >> T.C.y;
    cout<<'\n';
    if(areVerticesIntegerPoints(T)==true)
    {
    cout<<tr("Wierzcho³ki trójk¹ta s¹ punktami kratowymi")<<'\n';
    cout<<tr("Wewn¹trz trójk¹ta znajduje siê ")<<getInternalCount(T)<<tr(" punktów kratowych")<<'\n';
    cout<<tr("Na obwodzie trójk¹ta znajduje siê ")<<integersPointsOnSides(T)<<tr(" punktów kratowych")<<'\n';
    }
    line AB,BC,AC; /// Deklaruje istnienie prostych AB,BC i CA
    AB.P1 = T.A;
    AB.P2 = T.B;
    BC.P1 = T.B;
    BC.P2 = T.C;
    AC.P1 = T.A;
    AC.P2 = T.C;
   if(!Area(T))
    {
        bool p = operator==(T.A,T.B);
        bool q = operator==(T.A,T.C);
        bool r = operator==(T.B,T.C);
        if(p || q || r == true)
        {
            cout<<tr("Co najmniej dwa wprowadzone punkty siê pokrywaj¹ a zatem nie wyznaczaj¹ trójk¹ta!")<<'\n';
            cout<<tr("Proszê ponownie wprowadziæ prawid³owe wspó³rzêdne! ")<<'\n';
            cout<<'\n';
        }
        else{
        cout<<tr("Wprowadzone punkty nie wyznaczaj¹ trójk¹ta - s¹ wspó³liniowe!")<<'\n';
        cout<<tr("Proszê ponownie wprowadziæ prawid³owe wspó³rzêdne! ")<<'\n';
        cout<<'\n';
        }

    }
   else{
    Point O = CircleThroughThreePoints(T);
    Point OA = wektor(O,T.A);
    Point OB = wektor(O,T.B);
    Point OC = wektor(O,T.C);
    Point OH = Sum(OA,OB);
    OH = Sum(OH,OC);
    cout<<tr("Trójk¹t jest ")<<TypeOfTriangle(T)<<" i "<<TypeOfTriangle_Angles(T)<<'\n';
    cout<<tr("D³ugoœci boków wynosz¹: ")<<'\n';
    cout<<"AB = "<<LengthOfSegment(T.A,T.B)<<" BC = "<<LengthOfSegment(T.C,T.B)<<" AC = "<<LengthOfSegment(T.A,T.C)<<'\n';
    cout<<tr("Miary k¹tów (w stopniach) wynosz¹ odpowiednio:")<<" |<)BAC| = "<<(angels(T).p1)*(180.0/PI)<<" "<<" |<)ABC| = "<<(angels(T).p2)*(180.0/PI)<<" "<<" |<)BCA| = "<<(angels(T).p3)*(180.0/PI)<<""<<'\n';
    cout<<tr("Œrodki boków maj¹ wspó³rzêdne: ")<<'\n';
    cout<<"M_AB = ("<<MidP(T.A,T.B).x <<","<<MidP(T.A,T.B).y<<")"<<'\n';
    cout<<"M_BC = ("<<MidP(T.B,T.C).x<<","<<MidP(T.B,T.C).y<<")"<<'\n';
    cout<<"M_AC = ("<<MidP(T.A,T.C).x<<","<<MidP(T.A,T.C).y<<")"<<'\n';
    cout<<'\n';
    cout<<tr("Pole trójkata wynosi S = ")<<Area(T)<<'\n';
    cout<<tr("Obwód trójk¹ta wynosi L = ")<<Perimeter(T)<<'\n';
    cout<<'\n';
    cout<<tr("Promieñ okrêgu wpisanego w trójk¹t wynosi r = ")<<RadiusOfInCircle(T)<<'\n';
    cout<<tr("Promieñ okrêgu opisanego na trójk¹cie wynosi R = ")<<RadiusOfCircumCircle(T)<<'\n';
    cout<<'\n';
    cout<<tr("D³ugoœci wysokoœci trójk¹ta wynosz¹ odpowiednio: ")<<"h_a = "<<DistanceFromLine(BC,T.A)<<" h_b = "<<DistanceFromLine(AC,T.B)<<" h_c = "<<DistanceFromLine(AB,T.C)<<'\n';
    cout<<tr("D³ugoœci œrodkowych trójk¹ta wynosz¹ odpowiednio: ")<<"m_a = "<<LengthOfSegment(T.A,MidP(T.B,T.C))<<" m_b = "<<LengthOfSegment(T.B,MidP(T.A,T.C))<<" m_c = "<<LengthOfSegment(T.C,MidP(T.B,T.A))<<'\n';
    cout<<tr("D³ugoœci dwusiecznych trójk¹ta wynosz¹ odpowiednio: ")<<"d_a = "<<LengthOfAngelBisector(T).p1<<" d_b = "<<LengthOfAngelBisector(T).p2<<" d_c = "<<LengthOfAngelBisector(T).p3<<'\n';
    cout<<'\n';
    cout<<"Punkt G = ("<<Centroid(T).x<<","<<Centroid(T).y<<tr(") jest œrodkiem ciê¿koœci trójk¹ta - przeciêciem jego œrodkowych")<<'\n';
    cout<<"Punkt O = ("<<O.x<<","<<O.y<<tr(") jest œrodkiem okrêgu opisanego na trójk¹cie - przeciêciem symetralnych jego boków")<<'\n';

    cout<<"Punkt H = ("<<Orthocenter(T).x<<","<<Orthocenter(T).y<<tr(") jest ortocentrum trójk¹ta - przeciêciem jego wysokoœci")<<'\n';
    cout<<"Punkt I = ("<<Incenter(T).x<<","<<Incenter(T).y<<tr(") jest œrodkiem okrêgu wpisanego w trójk¹t - przeciêciem jego dwusiecznych")<<'\n';
    cout<<"Punkt L = ("<<Lemoine_Point(T).x<<","<<Lemoine_Point(T).y<<tr(") jest punktem Lemoine'a trójk¹ta - przeciêciem jego symedian")<<'\n';
    cout<<"Punkt G_e = ("<<Gergonne_Point(T).x<<","<<Gergonne_Point(T).y<<tr(") jest punktem Gergonne'a trójk¹ta")<<'\n';
    cout<<"Punkt N_a = ("<<Nagel_Point(T).x<<","<<Nagel_Point(T).y<<tr(") jest punktem Nagela'a trójk¹ta")<<'\n';


    }
    getchar();


    return 0;
}
 /// Bug nr 1 - czasem wspolrzedne O liczone przez funkcje FindCircumCenter s¹ niepoprawne (nie zgadza sie znak)
