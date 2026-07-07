#include <vector>
#include <variant>
#include <unordered_map>
#include <string>
#include <cmath>

namespace geo {

    double GEO_RADIUS = 6371000.0;
    double PI = std::acos(-1.0);

    struct GeoPoint {
        /// @brief WGS84 координаты
        double lat, lon, alt = 0.0;
    };

    class GeoLineString : public std::vector<GeoPoint> {
        /// @brief Гео-линия из N точек
    public:
        double length() const {
            /// @brief Вычисление длины линии с помощью формулы Гаверсинусов
            double dist = 0.0;

            for (size_t i = 1; i < this->size(); ++i) {
                double lat1 = this->at(i - 1).lat * PI / 180.0;
                double lon1 = this->at(i - 1).lon * PI / 180.0;
                double lat2 = this->at(i).lat * PI / 180.0;
                double lon2 = this->at(i).lon * PI / 180.0;
                double dlat = lat2 - lat1;
                double dlon = lon2 - lon1;
                double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
                        std::cos(lat1) * std::cos(lat2) *
                        std::sin(dlon / 2) * std::sin(dlon / 2);
                double c = 2 * std::atan2(std::sqrt(a), std::sqrt(1 - a));
                dist = dist + GEO_RADIUS * c;
            }
            return dist;
        };
    };

    class GeoPolygon {
        /// @brief Гео-полигон, фигура с внешним контуром и внутренними непересекающимися контурами
    public:
        bool contains(const GeoPoint& p) {
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
        bool pointInPolygon(const GeoPoint& p, GeoLineString line) {
            /// @brief Алгоритм ray-casting проверяет находится ли точка внутри замкнутой линии
            bool inside = false;

            for (size_t i = 0, j = line.size() - 1; i < line.size(); j = i++) {
                if (((line[i].lat > p.lat) != (line[j].lat > p.lat)) &&
                    (p.lon < (line[j].lon - line[i].lon) * (p.lat - line[i].lat) / 
                    (line[j].lat - line[i].lat) + line[i].lon)) {
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