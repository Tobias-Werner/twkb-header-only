//
// Created by Tobias Werner
//

#ifndef TWKB_HEADER_ONLY_TWKB_H
#define TWKB_HEADER_ONLY_TWKB_H

#include <vector>

namespace TWKB {

    using namespace std;

    using bytes_t = vector<unsigned char>;

    namespace GeometryType {

        enum Type {
            Point = 0x01,
            Line = 0x02,
            Polygon = 0x03,
            MultiPoint = 0x04,
            MultiLineString = 0x05,
            MultiPolygon = 0x06,
            GeometryCollection = 0x07
        };
    }

    class Geometry {

        // Structure
        //
        // type_and_prec     byte
        //      Bit | Description
        //      ----|------------
        //      1-4 | Geometry type
        //      5-8 | Geometry Precision
        // metadata_header   byte
        // [extended_dims]   byte
        // [size]            uvarint
        // [bounds]          bbox
        // pointarray        varint[]



    public:
        bytes_t getBytes() {
            return bytes;
        }

        Geometry(GeometryType::Type type, unsigned char precision) : bytes({0x00, 0x00}) {
            bytes[0] |= encodeZigZag(precision);
            bytes[0] <<= 4;
            bytes[0] |= type;
        }

        bytes_t bytes;

        static unsigned int decodeVarint(bytes_t bytes) {

            unsigned int result = 0x00;

            for (size_t i = bytes.size(); i > 0; i--) {
                result |= (bytes.at(i - 1) & 0x7F);
                result <<= (i - 1 > 0) ? 7 : 0;
            }

            return result;
        }

        static bytes_t encodeVarint(unsigned int value) {
            bytes_t result;

            do {
                result.push_back(0x80 | (value & 0x7F));
                value >>= 7;
            } while (value > 0x00);

            result.back() &= 0x7F;

            return result;
        }

        static unsigned int encodeZigZag(int value) {
            return (value << 1) ^ (value >> 31);
        }

        static int decodeZigZag(unsigned int value) {
            return -(value & 1) ^ (value >> 1);
        }

    };

    class Point : public Geometry {
    public:
        Point(double lat, double lon, unsigned char precision) : Geometry(GeometryType::Type::Point, precision) {

        }
    };

    class Line : public Geometry {
    public:
        Line(double lat, double lon, unsigned char precision) : Geometry(GeometryType::Type::Point, precision) {

        }
    };

    class Polygon : public Geometry {
    public:
        Polygon(double lat, double lon, unsigned char precision) : Geometry(GeometryType::Type::Point, precision) {

        }
    };


}

#endif //TWKB_HEADER_ONLY_TWKB_H
