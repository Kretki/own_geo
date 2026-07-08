#include <vector>
#include <variant>
#include <unordered_map>
#include <string>
#include <cmath>
#include <stdexcept>
#include <numbers>

namespace geo {

    inline constexpr double EARTH_RADIUS = 6371007.2;
    inline constexpr double PI = std::numbers::pi;
    inline constexpr double a  = 6378137.0;
    inline constexpr double f  = 1.0 / 298.257223563;
    inline constexpr double b  = a * (1.0 - f);
    inline constexpr double e2 = f * (2.0 - f);

    class GeoPoint {
        /// @brief WGS84 и ECEF координаты
    public:
        GeoPoint() : lat(0.0), lon(0.0), alt(0.0) {};
        static GeoPoint fromWGS(double latitude, double longitude, double altitude) {
            auto p = GeoPoint();
            p.setWGS(latitude, longitude, altitude);
            return p;
        };
        static GeoPoint fromECEF(double X, double Y, double Z) {
            auto p = GeoPoint();
            p.setECEF(X, Y, Z);
            return p;
        }

        GeoPoint operator-(const GeoPoint& other) const {
            return GeoPoint::fromECEF(this->X - other.X, this->Y - other.Y, this->Z - other.Z);
        }
        bool operator==(const GeoPoint& other) const {
            constexpr double latLonEps = 1e-6;
            constexpr double altEps = 1e-4;
            return (std::abs(lat - other.lat) < latLonEps) &&
                (std::abs(lon - other.lon) < latLonEps) &&
                (std::abs(alt - other.alt) < altEps);
        }
        bool operator!=(const GeoPoint& other) const {
            return !(*this == other);
        }
        
        void setLat(double latitude) { 
            if (std::abs(latitude) > 90.0) throw std::invalid_argument("Latitude out of range");
            lat = latitude;
            this->fillECEF();
        }
        void setLon(double longitude) { 
            if (std::abs(longitude) > 180.0) throw std::invalid_argument("Longitude out of range");
            lon = longitude;
            this->fillECEF();
        }
        void setAlt(double altitude) {
            alt = altitude;
            this->fillECEF();
        }
        void setWGS(double latitude, double longitude, double altitude) {
            if (std::abs(latitude) > 90.0) throw std::invalid_argument("Latitude out of range");
            if (std::abs(longitude) > 180.0) throw std::invalid_argument("Longitude out of range");
            lon = longitude;
            lat = latitude;
            alt = altitude;
            this->fillECEF();
        }

        void setX(double X) { this->X = X; fillWGS(); }
        void setY(double Y) { this->Y = Y; fillWGS(); }
        void setZ(double Z) { this->Z = Z; fillWGS(); }
        void setECEF(double X, double Y, double Z) {
            this->X = X;
            this->Y = Y;
            this->Z = Z;
            this->fillWGS();
        }

        double getLat() const { return lat; }
        double getLon() const { return lon; }
        double getAlt() const { return alt; }
        double getX() const { return X; }
        double getY() const { return Y; }
        double getZ() const { return Z; }

        double distanceTo(const GeoPoint& other) const {
            double horizontal = haversineDistance(other);
            double vertical = other.alt - alt;
            return std::sqrt(horizontal * horizontal + vertical * vertical);
        }
    private:
        //WGS84 coords
        double lat, lon, alt; 
        //Geocentric coords
        double X, Y, Z;

        static double toRadians(double degrees) {
            return degrees * PI / 180.0;
        }
        static double toDegrees(double radians) {
            return radians * 180.0 / PI;
        }

        double fillECEF() {
            const double sinLat = std::sin(toRadians(this->lat));
            const double cosLat = std::cos(toRadians(this->lat));
            const double sinLon = std::sin(toRadians(this->lon));
            const double cosLon = std::cos(toRadians(this->lon));
            const double N = a / std::sqrt(1.0 - e2 * sinLat * sinLat);
            this->X = (N + this->alt) * cosLat * cosLon;
            this->Y = (N + this->alt) * cosLat * sinLon;
            this->Z = (N * (1.0 - e2) + this->alt) * sinLat;
        }
        double fillWGS() {
            this->lon = toDegrees(std::atan2(this->Y, this->X));
            const double p = std::sqrt(this->X * this->X + this->Y * this->Y);
            if (p < 1e-9) {
                this->lat = toDegrees((this->Z >= 0.0) ? (std::acos(-1.0) * 0.5)
                                                       : (std::acos(-1.0) * 0.5));
                this->alt = std::fabs(this->Z) - b;
            }
        }

        double haversineDistance(const GeoPoint& other) const {
            /// @brief Вычисление расстояния по шару с помощью формулы Гаверсинусов
            //TODO Переделать на WGS84
            double lat1 = toRadians(lat);
            double lon1 = toRadians(lon);
            double lat2 = toRadians(other.lat);
            double lon2 = toRadians(other.lon);

            double dLat = lat2 - lat1;
            double dLon = lon2 - lon1;

            double a = std::sin(dLat / 2) * std::sin(dLat / 2) +
                    std::cos(lat1) * std::cos(lat2) *
                    std::sin(dLon / 2) * std::sin(dLon / 2);

            double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
            return EARTH_RADIUS * c;
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
        bool pointInPolygon(const GeoPoint& p, GeoLineString line) const {
            /// @brief Алгоритм ray-casting проверяет находится ли точка внутри замкнутой линии
            //TODO Переделать учет переход через -+180 градусов
            bool inside = false;
            for (size_t i = 0, j = line.size() - 1; i < line.size(); j = i++) {
                if (((line[i].getLat() > p.getLat()) != (line[j].getLat() > p.getLat())) &&
                    (p.getLon() < (line[j].getLon() - line[i].getLon()) * 
                    (p.getLat() - line[i].getLat()) / 
                    (line[j].getLat() - line[i].getLat()) + 
                    line[i].getLon())) {
                        inside = !inside;
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