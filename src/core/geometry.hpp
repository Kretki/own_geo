#include <vector>
#include <variant>
#include <unordered_map>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numbers>

namespace geo {
    inline constexpr double PI = std::numbers::pi;
    inline constexpr double a  = 6378137.0;
    inline constexpr double f  = 1.0 / 298.257223563;
    inline constexpr double b  = a * (1.0 - f);
    inline constexpr double e2 = f * (2.0 - f);
    inline constexpr double ecefEps = 1e-4;
    inline constexpr double wgsEps = 1e-12;
    inline constexpr double coincidentEps = 1e-24;

    class GeoPoint {
        /// @brief WGS84 и ECEF координаты
        /// @concept Одновременное вычисление ECEF и WGS84 координат позволяет избежать больших пересчетов
    public:
        GeoPoint() {
            /// @brief Основной конструктор - по-умолчанию
            setWGS(0.0, 0.0, 0.0); 
        };
        static GeoPoint fromWGS(double latitude, double longitude, double altitude) {
            /// @brief Создание по WGS84 с автоматическим заполнением ECEF
            auto p = GeoPoint();
            p.setWGS(latitude, longitude, altitude);
            return p;
        };
        static GeoPoint fromECEF(double X, double Y, double Z) {
            /// @brief Создание по ECEF с автоматическим заполнением WGS84
            auto p = GeoPoint();
            p.setECEF(X, Y, Z);
            return p;
        }

        GeoPoint operator-(const GeoPoint& other) const {
            /// @brief Разность векторов в ECEF
            return GeoPoint::fromECEF(_X - other._X, _Y - other._Y, _Z - other._Z);
        }
        bool operator==(const GeoPoint& other) const {
            /// @brief Сравнение точек
            /// @important помним про то, что double может не равняться точно
            return (std::abs(_X - other._X) < ecefEps) &&
                   (std::abs(_Y - other._Y) < ecefEps) &&
                   (std::abs(_Z - other._Z) < ecefEps);
        }
        bool operator!=(const GeoPoint& other) const {
            /// @brief Сравнение точек
            return !(*this == other);
        }
        
        void setLat(double latitude) { 
            if (std::abs(latitude) > 90.0) throw std::invalid_argument("Latitude out of range");
            _lat = latitude;
            fillECEF();
        }
        void setLon(double longitude) { 
            if (std::abs(longitude) > 180.0) throw std::invalid_argument("Longitude out of range");
            _lon = longitude;
            fillECEF();
        }
        void setAlt(double altitude) {
            _alt = altitude;
            fillECEF();
        }
        void setWGS(double latitude, double longitude, double altitude) {
            if (std::abs(latitude) > 90.0) throw std::invalid_argument("Latitude out of range");
            if (std::abs(longitude) > 180.0) throw std::invalid_argument("Longitude out of range");
            _lon = longitude;
            _lat = latitude;
            _alt = altitude;
            //После изменения всех WGS заполняем ECEF
            fillECEF();
        }

        void setX(double X) { _X = X; fillWGS(); }
        void setY(double Y) { _Y = Y; fillWGS(); }
        void setZ(double Z) { _Z = Z; fillWGS(); }
        void setECEF(double X, double Y, double Z) {
            _X = X;
            _Y = Y;
            _Z = Z;
            //После изменения всех ECEF заполняем WGS
            fillWGS();
        }

        double getLat() const { return _lat; }
        double getLon() const { return _lon; }
        double getAlt() const { return _alt; }
        double getX() const { return _X; }
        double getY() const { return _Y; }
        double getZ() const { return _Z; }

        double distanceTo(const GeoPoint& other) const {
            /// @brief Расстояние по WGS84
            double horizontal = vincentyDistance(other);
            double vertical = other._alt - _alt;
            return std::sqrt(horizontal * horizontal + vertical * vertical);
        }
    private:
        //WGS84 coords
        double _lat, _lon, _alt; 
        //Geocentric coords
        double _X, _Y, _Z;

        static double toRadians(double degrees) {
            return degrees * PI / 180.0;
        }
        static double toDegrees(double radians) {
            return radians * 180.0 / PI;
        }

        void fillECEF() {
            /// @brief Заполнение ECEF после изменения WGS84
            const double sinLat = std::sin(toRadians(_lat));
            const double cosLat = std::cos(toRadians(_lat));
            const double sinLon = std::sin(toRadians(_lon));
            const double cosLon = std::cos(toRadians(_lon));
            const double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
            _X = (N + _alt) * cosLat * cosLon;
            _Y = (N + _alt) * cosLat * sinLon;
            _Z = (N * (1.0 - e2) + _alt) * sinLat;
        }
        void fillWGS() {
            /// @brief Заполнение WGS84 после изменения ECEF
            _lon = toDegrees(std::atan2(_Y, _X));
            const double p = std::sqrt(_X * _X + _Y * _Y);
            if (p < wgsEps) {
                _lat = toDegrees((_Z >= 0.0) ? 90.0 : -90.0);
                _alt = std::fabs(_Z) - b;
            }
            else {
                _lat = std::atan2(_Z, p * (1.0 - e2));
                _alt = 0.0;
                for (int i = 0; i < 10; ++i) {
                    const double sinLat = std::sin(toRadians(_lat));
                    const double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
                    _alt = p / std::cos(_lat) - N;
                    const double newLat = std::atan((_Z / p) / (1.0 - e2 * N / (N + _alt)));
                    if (std::abs(newLat - _lat) < wgsEps) break;
                    _lat = newLat;
                }
            }
        }

        double vincentyDistance(const GeoPoint& other) const {
            /// @brief Вычисление расстояния по шару с помощью формулы Винценти по WGS84
            double lat1 = toRadians(_lat);
            double lon1 = toRadians(_lon);
            double lat2 = toRadians(other._lat);
            double lon2 = toRadians(other._lon);

            double L = lon2 - lon1;

            double tanU1 = (1.0 - f) * std::tan(lat1);
            double cosU1 = 1.0 / std::sqrt(1.0 + tanU1 * tanU1);
            double sinU1 = tanU1 * cosU1;

            double tanU2 = (1.0 - f) * std::tan(lat2);
            double cosU2 = 1.0 / std::sqrt(1.0 + tanU2 * tanU2);
            double sinU2 = tanU2 * cosU2;

            double lambda = L;
            double sinLambda, cosLambda, sinSigma, cosSigma, sigma, cos2SigmaM, cosSqAlpha, lambdaP;

            do {
                sinLambda = std::sin(lambda);
                cosLambda = std::cos(lambda);

                double sinSqSigma = (cosU2 * sinLambda) * (cosU2 * sinLambda) +
                                    (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) *
                                    (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda);

                if (std::abs(sinSqSigma) < coincidentEps) {
                    return 0.0;
                }

                sinSigma = std::sqrt(sinSqSigma);
                cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
                sigma = std::atan2(sinSigma, cosSigma);

                double sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
                cosSqAlpha = 1.0 - sinAlpha * sinAlpha;

                cos2SigmaM = (cosSqAlpha != 0.0)
                    ? (cosSigma - 2.0 * sinU1 * sinU2 / cosSqAlpha)
                    : 0.0;

                double C = f / 16.0 * cosSqAlpha * (4.0 + f * (4.0 - 3.0 * cosSqAlpha));

                lambdaP = lambda;
                lambda = L + (1.0 - C) * f * sinAlpha *
                        (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));
            } while (std::abs(lambda - lambdaP) > wgsEps);

            if (std::abs(lambda) > PI) {
                return std::numeric_limits<double>::quiet_NaN();
            }

            double uSq = cosSqAlpha * (a * a - b * b) / (b * b);
            double A = 1.0 + uSq / 16384.0 * (4096.0 + uSq * (-768.0 + uSq * (320.0 - 175.0 * uSq)));
            double B = uSq / 1024.0 * (256.0 + uSq * (-128.0 + uSq * (74.0 - 47.0 * uSq)));

            double deltaSigma = B * sinSigma * (
                cos2SigmaM + B / 4.0 * (
                    cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM) -
                    B / 6.0 * cos2SigmaM * (-3.0 + 4.0 * sinSigma * sinSigma) *
                    (-3.0 + 4.0 * cos2SigmaM * cos2SigmaM)
                )
            );

            double s = b * A * (sigma - deltaSigma);
            return s;
        }
    };

    class GeoLineString {
        /// @brief Гео-линия из N точек
    public:
        explicit GeoLineString(std::vector<GeoPoint> pts) : points(std::move(pts)) {};
        const GeoPoint& operator[](size_t index) const {
            if (index >= points.size()) throw std::out_of_range("Index out of range");
            return points[index];
        }
        double length() const {
            double dist = 0.0;

            for (size_t i = 1; i < points.size(); ++i) {
                dist = dist + points[i].distanceTo(points[i - 1]);
            }
            return dist;
        };
        int size() const {
            return points.size();
        }
        void push_back(const GeoPoint& p) {
            points.push_back(p);
        }
    private:
        std::vector<GeoPoint> points;
    };

    class GeoPolygon {
        /// @brief Гео-полигон, фигура с внешним контуром и внутренними непересекающимися контурами
    public:
        bool contains(const GeoPoint& p) const {
            if (!isValid()) return false;
            if (!pointInPolygon(p, outerRing)) return false;
            for (auto innerRing : innerRings) {
                if (pointInPolygon(p, innerRing)) return false;
            }
            return true;
        };
        bool isValid() const {
            if (outerRing.size() < 3) return false;
            for (auto innerRing : innerRings) {
                if (innerRing.size() < 3) return false;
            }
            return true;
        };
    private:
        bool pointInPolygon(const GeoPoint& p, const GeoLineString& line) const {
            /// @brief Алгоритм ray-casting проверяет находится ли точка внутри замкнутой линии
            bool inside = false;
            auto normLon = [](double lon) -> double {
                return (lon < 0.0) ? lon + 360.0 : lon;
            };
            double pLat = p.getLat();
            double pLonNorm = normLon(p.getLon());
            for (size_t i = 0, j = line.size() - 1; i < line.size(); j = i++) {
                double lat_i = line[i].getLat();
                double lat_j = line[j].getLat();
                if (((lat_i > pLat) != (lat_j > pLat))) {
                    double lon_i = normLon(line[i].getLon());
                    double lon_j = normLon(line[j].getLon());
                    double x_intersect;
                    double denom = lat_j - lat_i;
                    if (std::abs(denom) < wgsEps) {
                        x_intersect = lon_i;
                    } else {
                        x_intersect = lon_i + (lon_j - lon_i) *
                                    (pLat - lat_i) / denom;
                    }
                    if (pLonNorm < x_intersect) {
                        inside = !inside;
                    }
                }
            }
            return inside;
        }

        GeoLineString outerRing;
        std::vector<GeoLineString> innerRings;
    };

    class GeoFeature {
        /// @brief Инкапсулирует гео-данные и мета-данные
        std::variant<GeoPoint, GeoLineString, GeoPolygon> geometry;
        std::unordered_map<std::string, std::variant<std::string, double, int>> properties;
    };
}