//
// Created by Tobias Werner
//

#ifndef TWKB_HEADER_ONLY_TWKB_H
#define TWKB_HEADER_ONLY_TWKB_H

#include <list>

namespace TWKB {

    using namespace std;

    using bytes_t = list<u_char>;

    struct LocationXY {
        double x, y;

        LocationXY(double x, double y) : x(x), y(y) {}
    };

    struct LocationXYZ {
        double x, y, z;

        LocationXYZ(double x, double y, double z) : x(x), y(y), z(z) {}
    };

    struct LocationXYZT {
        double x, y, z, t;

        LocationXYZT(double x, double y, double z, double t) : x(x), y(y), z(z), t(t) {}
    };


    class GeomFactory {

    private:

        static int32_t shrink(const double &value, const signed char &precision) {

            const bool negative = (precision < 0);
            u_char p = (negative) ? static_cast<u_char>(-precision) : static_cast<u_char>(precision);
            auto unrounded = (negative) ? pow10(p) / value : pow10(p) * value;
            return static_cast<int32_t>(round(unrounded));

        }

        static bytes_t encode(const double &value, const signed char &precision) {
            return encodeVarint(encodeZigZag(shrink(value, precision)));
        }

        static void append(bytes_t &twkb) {

        }


        template<typename T, typename... Types>
        static void append(bytes_t &twkb, const T &var1, const Types... var2) {
            twkb.insert(twkb.end(), var1.begin(), var1.end());
            append(twkb, var2...);
        }


        static uint32_t pow10(u_char exponent) {
            static uint32_t pow10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
            return pow10[exponent];
        }

        static int32_t decodeVarint(const bytes_t &bytes) {

            int32_t result = 0x00;

            size_t i = bytes.size();
            for (auto elem = bytes.rbegin(); elem != bytes.rend(); ++elem, i--) {
                result |= (*elem & 0x7F);
                result <<= (i - 1 > 0) ? 7 : 0;
            }

            return result;
        }

        static bytes_t encodeVarint(uint32_t value) {
            bytes_t result;

            do {
                result.push_back(0x80 | (value & 0x7F));
                value >>= 7;
            } while (value > 0x00);

            result.back() &= 0x7F;

            return result;
        }

        static uint32_t encodeZigZag(const int32_t &value) {
            return static_cast<uint32_t >((value << 1) ^ (value >> 31));
        }

        static int32_t decodeZigZag(const uint32_t &value) {
            return -(value & 1) ^ (value >> 1);
        }

        static bytes_t createHeader(const u_char &type, const signed char &precisionXY, const bool &bbox, const bool &sizes, const bool &idlist, const bool &extendedDimensions,
                                    const bool &emptyGeom) {

            u_char header = type | encodeZigZag(precisionXY) << 4;
            bytes_t twkb({header});

            u_char metadataHeader = 0x00;
            metadataHeader |= (bbox) ? 0x01 : 0x00;
            metadataHeader |= (sizes) ? 0x02 : 0x00;
            metadataHeader |= (idlist) ? 0x04 : 0x00;
            metadataHeader |= (extendedDimensions) ? 0x08 : 0x00;
            metadataHeader |= (emptyGeom) ? 0x10 : 0x00;
            twkb.push_back(metadataHeader);

            return twkb;
        }

        static bytes_t createHeader(const u_char &type, const signed char &precisionXY, const bool &bbox, const bool &sizes, const bool &idlist, const bool &extendedDimensions,
                                    const bool &emptyGeom, const u_char &precisionZ) {

            auto twkb = createHeader(type, precisionXY, bbox, sizes, idlist, extendedDimensions, emptyGeom);
            u_char extendedInformation = 0x01;
            extendedInformation |= (precisionZ << 2);
            twkb.push_back(extendedInformation);
            return twkb;
        }

        static bytes_t createHeader(const u_char &type, const signed char &precisionXY, const bool &bbox, const bool &sizes, const bool &idlist, const bool &extendedDimensions,
                                    const bool &emptyGeom, const u_char &precisionZ, const u_char &precisionT) {

            auto twkb = createHeader(type, precisionXY, bbox, sizes, idlist, extendedDimensions, emptyGeom);
            u_char extendedInformation = 0x03;
            extendedInformation |= (precisionZ << 2);
            extendedInformation |= (precisionT << 5);
            twkb.push_back(extendedInformation);
            return twkb;
        }

    public:

        bytes_t makePoint(const LocationXY &location, const signed char &precisionXY) {

            bytes_t twkb = createHeader(0x01, precisionXY, false, false, false, false, false);

            // Locations
            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);

            append(twkb, x, y);

            return twkb;
        }

        bytes_t makePoint(const LocationXYZ &location, const signed char &precisionXY, const u_char &precisionZ) {

            bytes_t twkb = createHeader(0x01, precisionXY, false, false, false, true, false, precisionZ);

            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);
            auto z = encode(location.z, static_cast<signed char>(precisionZ));

            append(twkb, x, y, z);

            return twkb;
        }

        bytes_t makePoint(const LocationXYZT &location, const signed char &precisionXY, const u_char &precisionZ, const u_char &precisionT) {

            bytes_t twkb = createHeader(0x01, precisionXY, false, false, false, true, false, precisionZ, precisionT);

            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);
            auto z = encode(location.z, static_cast<signed char>(precisionZ));
            auto t = encode(location.t, static_cast<signed char>(precisionT));

            append(twkb, x, y, z, t);

            return twkb;
        }

        bytes_t makeLine(const vector<LocationXY> &locations, const signed char &precisionXY, const bool &bbox) {

            bytes_t twkb = createHeader(0x02, precisionXY, bbox, false, false, false, false);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;

                for (const auto &location : locations) {
                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));

                append(twkb, bboxX, deltaX, bboxY, deltaY);

            }

            // Set number of containing points
            bytes_t npoints = encodeVarint(static_cast<uint32_t >(locations.size()));

            append(twkb, npoints);

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));

            append(twkb, x, y);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));

                append(twkb, x, y);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);

            }

            return twkb;
        }

        bytes_t makeLine(const vector<LocationXYZ> &locations, const signed char &precisionXY, const u_char &precisionZ, const bool &bbox) {

            bytes_t twkb = createHeader(0x02, precisionXY, bbox, false, false, true, false, precisionZ);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;

                for (const auto &location : locations) {
                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpZ < minz) minz = tmpZ;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                    if (tmpZ > maxz) maxz = tmpZ;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ);

            }


            // Set number of containing points
            bytes_t npoints = encodeVarint(static_cast<uint32_t >(locations.size()));
            twkb.insert(twkb.end(), npoints.begin(), npoints.end());

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, static_cast<signed char>(precisionZ));

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));

            append(twkb, x, y, z);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, static_cast<signed char>(precisionZ)) - zShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));

                append(twkb, x, y, z);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, static_cast<signed char>(precisionZ));

            }

            return twkb;
        }

        bytes_t
        makeLine(const vector<LocationXYZT> &locations, const signed char &precisionXY, const u_char &precisionZ, const u_char &precisionT, const bool &bbox) {

            bytes_t twkb = createHeader(0x02, precisionXY, bbox, false, false, true, false, precisionZ, precisionT);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;
                int32_t mint = shrink(locations.front().t, static_cast<signed char>(precisionT));
                int32_t maxt = mint;

                for (const auto &location : locations) {
                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));
                    int32_t tmpT = shrink(location.t, static_cast<signed char>(precisionT));

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpZ < minz) minz = tmpZ;
                    if (tmpT < mint) mint = tmpT;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                    if (tmpZ > maxz) maxz = tmpZ;
                    if (tmpT > maxt) maxt = tmpT;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t bboxT = encodeVarint(encodeZigZag(mint));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));
                bytes_t deltaT = encodeVarint(encodeZigZag(maxt - mint));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ, bboxT, deltaT);

            }


            // Set number of containing points
            bytes_t npoints = encodeVarint(static_cast<uint32_t >(locations.size()));
            twkb.insert(twkb.end(), npoints.begin(), npoints.end());

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, static_cast<signed char>(precisionZ));
            int32_t tShrinked = shrink(locations.front().t, static_cast<signed char>(precisionT));

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));
            bytes_t t = encodeVarint(encodeZigZag(tShrinked));

            append(twkb, x, y, z, t);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, static_cast<signed char>(precisionZ)) - zShrinked;
                int32_t deltaT = shrink(locations[i].t, static_cast<signed char>(precisionT)) - tShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));
                t = encodeVarint(encodeZigZag(deltaT));

                append(twkb, x, y, z, t);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, static_cast<signed char>(precisionZ));
                tShrinked = shrink(locations[i].t, static_cast<signed char>(precisionT));

            }

            return twkb;
        }

        bytes_t makePolygon(const vector<vector<LocationXY>> &rings, const signed char &precisionXY, const bool &bbox) {

            bytes_t twkb = createHeader(0x03, precisionXY, bbox, false, false, false, false);

            if (bbox) {

                int32_t minx = shrink(rings.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(rings.front().front().y, precisionXY);
                int32_t maxy = miny;

                for (const auto &ring : rings) {

                    for (const auto &location : ring) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);

                        if (tmpX < minx) minx = tmpX;
                        if (tmpY < miny) miny = tmpY;
                        if (tmpX > maxx) maxx = tmpX;
                        if (tmpY > maxy) maxy = tmpY;
                    }

                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));

                append(twkb, bboxX, deltaX, bboxY, deltaY);

            }


            // Insert number of containing rings
            bytes_t nrings = encodeVarint(static_cast<uint32_t >(rings.size()));
            append(twkb, nrings);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;

            for (const auto &pointsInRing: rings) {

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(static_cast<uint32_t >(pointsInRing.size()));
                append(twkb, npoints);

                for (const auto &point : pointsInRing) {

                    int32_t currentXFull = shrink(point.x, precisionXY);
                    int32_t currentYFull = shrink(point.y, precisionXY);

                    int32_t deltaX = currentXFull - lastXFull;
                    int32_t deltaY = currentYFull - lastYFull;

                    lastXFull = currentXFull;
                    lastYFull = currentYFull;

                    bytes_t deltaXAsZigZag = encodeVarint(encodeZigZag(deltaX));
                    bytes_t deltaYAsZigZag = encodeVarint(encodeZigZag(deltaY));

                    append(twkb, deltaXAsZigZag, deltaYAsZigZag);

                }
            }


            return
                    twkb;
        }

        bytes_t

        makePolygon(const vector<vector<LocationXYZ>> &rings, const signed char &precisionXY, const u_char &precisionZ, const bool &bbox) {

            bytes_t twkb = createHeader(0x03, precisionXY, bbox, false, false, true, false, precisionZ);

            if (bbox) {

                int32_t minx = shrink(rings.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(rings.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(rings.front().front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;

                for (const auto &ring : rings) {

                    for (const auto &location : ring) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));

                        if (tmpX < minx) minx = tmpX;
                        if (tmpY < miny) miny = tmpY;
                        if (tmpZ < minz) minz = tmpZ;
                        if (tmpX > maxx) maxx = tmpX;
                        if (tmpY > maxy) maxy = tmpY;
                        if (tmpZ > maxz) maxz = tmpZ;
                    }

                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ);

            }

            // Insert number of containing rings
            bytes_t nrings = encodeVarint(static_cast<uint32_t>(rings.size()));
            append(twkb, nrings);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;

            for (const auto &pointsInRing : rings) {

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(static_cast<uint32_t >(pointsInRing.size()));
                append(twkb, npoints);

                for (const auto &point : pointsInRing) {

                    int32_t currentXFull = shrink(point.x, precisionXY);
                    int32_t currentYFull = shrink(point.y, precisionXY);
                    int32_t currentZFull = shrink(point.z, static_cast<signed char>(precisionZ));

                    int32_t deltaX = currentXFull - lastXFull;
                    int32_t deltaY = currentYFull - lastYFull;
                    int32_t deltaZ = currentZFull - lastZFull;

                    lastXFull = currentXFull;
                    lastYFull = currentYFull;
                    lastZFull = currentZFull;

                    bytes_t deltaXAsZigZag = encodeVarint(encodeZigZag(deltaX));
                    bytes_t deltaYAsZigZag = encodeVarint(encodeZigZag(deltaY));
                    bytes_t deltaZAsZigZag = encodeVarint(encodeZigZag(deltaZ));

                    append(twkb, deltaXAsZigZag, deltaYAsZigZag, deltaZAsZigZag);

                }
            }


            return twkb;
        }

        bytes_t
        makePolygon(const vector<vector<LocationXYZT>> &rings, const signed char &precisionXY, const u_char &precisionZ, const u_char &precisionT, const bool &bbox) {

            bytes_t twkb = createHeader(0x03, precisionXY, bbox, false, false, true, false, precisionZ, precisionT);

            if (bbox) {

                int32_t minx = shrink(rings.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(rings.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(rings.front().front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;
                int32_t mint = shrink(rings.front().front().t, static_cast<signed char>(precisionT));
                int32_t maxt = mint;

                for (const auto &ring : rings) {

                    for (const auto &location : ring) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));
                        int32_t tmpT = shrink(location.t, static_cast<signed char>(precisionT));

                        if (tmpX < minx) minx = tmpX;
                        if (tmpY < miny) miny = tmpY;
                        if (tmpZ < minz) minz = tmpZ;
                        if (tmpT < mint) mint = tmpT;
                        if (tmpX > maxx) maxx = tmpX;
                        if (tmpY > maxy) maxy = tmpY;
                        if (tmpZ > maxz) maxz = tmpZ;
                        if (tmpT > maxt) maxt = tmpT;
                    }

                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t bboxT = encodeVarint(encodeZigZag(mint));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));
                bytes_t deltaT = encodeVarint(encodeZigZag(maxt - mint));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ, bboxT, deltaT);

            }

            // Insert number of containing rings
            bytes_t nrings = encodeVarint(static_cast<uint32_t >(rings.size()));
            append(twkb, nrings);


            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;
            int32_t lastTFull = 0;

            for (const auto &pointsInRing : rings) {

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(static_cast<uint32_t >(pointsInRing.size()));
                append(twkb, npoints);

                for (const auto &point : pointsInRing) {

                    int32_t currentXFull = shrink(point.x, precisionXY);
                    int32_t currentYFull = shrink(point.y, precisionXY);
                    int32_t currentZFull = shrink(point.z, static_cast<signed char>(precisionZ));
                    int32_t currentTFull = shrink(point.t, static_cast<signed char>(precisionT));

                    int32_t deltaX = currentXFull - lastXFull;
                    int32_t deltaY = currentYFull - lastYFull;
                    int32_t deltaZ = currentZFull - lastZFull;
                    int32_t deltaT = currentTFull - lastTFull;

                    lastXFull = currentXFull;
                    lastYFull = currentYFull;
                    lastZFull = currentZFull;
                    lastTFull = currentTFull;

                    bytes_t deltaXAsZigZag = encodeVarint(encodeZigZag(deltaX));
                    bytes_t deltaYAsZigZag = encodeVarint(encodeZigZag(deltaY));
                    bytes_t deltaZAsZigZag = encodeVarint(encodeZigZag(deltaZ));
                    bytes_t deltaTAsZigZag = encodeVarint(encodeZigZag(deltaT));

                    append(twkb, deltaXAsZigZag, deltaYAsZigZag, deltaZAsZigZag, deltaTAsZigZag);

                }
            }


            return twkb;
        }

        bytes_t makeMultiPoint(const vector<LocationXY> &locations, const signed char &precisionXY, const bool &bbox) {

            bytes_t twkb = createHeader(0x04, precisionXY, bbox, false, false, false, false);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;

                for (const auto &location : locations) {

                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;

                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));

                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));

                append(twkb, bboxX, deltaX, bboxY, deltaY);

            }

            // Set number of containing points
            bytes_t npoints = encodeVarint(static_cast<uint32_t >(locations.size()));
            append(twkb, npoints);

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));

            append(twkb, x, y);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));

                append(twkb, x, y);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);

            }


            return twkb;
        }

        bytes_t makeMultiPoint(const vector<LocationXYZ> &locations, const signed char &precisionXY, const u_char &precisionZ, const bool &bbox) {

            bytes_t twkb = createHeader(0x04, precisionXY, bbox, false, false, true, false, precisionZ);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;

                for (const auto &location : locations) {

                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpZ < minz) minz = tmpZ;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                    if (tmpZ > maxz) maxz = tmpZ;

                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));

                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ);

            }

            // Set number of containing points
            bytes_t npoints = encodeVarint(static_cast<uint32_t >(locations.size()));
            append(twkb, npoints);

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, static_cast<signed char>(precisionZ));

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));

            append(twkb, x, y, z);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, static_cast<signed char>(precisionZ)) - zShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));

                append(twkb, x, y, z);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, static_cast<signed char>(precisionZ));

            }


            return twkb;
        }

        bytes_t
        makeMultiPoint(const vector<LocationXYZT> &locations, const signed char &precisionXY, const u_char &precisionZ, const u_char &precisionT, const bool &bbox) {

            bytes_t twkb = createHeader(0x04, precisionXY, bbox, false, false, true, false, precisionZ, precisionT);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;
                int32_t mint = shrink(locations.front().t, static_cast<signed char>(precisionT));
                int32_t maxt = mint;

                for (const auto &location : locations) {

                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));
                    int32_t tmpT = shrink(location.t, static_cast<signed char>(precisionT));

                    if (tmpX < minx) minx = tmpX;
                    if (tmpY < miny) miny = tmpY;
                    if (tmpZ < minz) minz = tmpZ;
                    if (tmpT < mint) mint = tmpT;
                    if (tmpX > maxx) maxx = tmpX;
                    if (tmpY > maxy) maxy = tmpY;
                    if (tmpZ > maxz) maxz = tmpZ;
                    if (tmpT > maxt) maxt = tmpT;
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t bboxT = encodeVarint(encodeZigZag(mint));

                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));
                bytes_t deltaT = encodeVarint(encodeZigZag(maxt - mint));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ, bboxT, deltaT);

            }

            // Set number of containing points
            bytes_t npoints = encodeVarint(static_cast<uint32_t >(locations.size()));
            append(twkb, npoints);

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, static_cast<signed char>(precisionZ));
            int32_t tShrinked = shrink(locations.front().t, static_cast<signed char>(precisionT));

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));
            bytes_t t = encodeVarint(encodeZigZag(tShrinked));

            append(twkb, x, y, z, t);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, static_cast<signed char>(precisionZ)) - zShrinked;
                int32_t deltaT = shrink(locations[i].t, static_cast<signed char>(precisionT)) - tShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));
                t = encodeVarint(encodeZigZag(deltaT));

                append(twkb, x, y, z, t);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, static_cast<signed char>(precisionZ));
                tShrinked = shrink(locations[i].t, static_cast<signed char>(precisionT));

            }


            return twkb;
        }

        bytes_t
        makeMultiLine(const vector<vector<LocationXYZT>> &lines, const signed char &precisionXY, const u_char &precisionZ, const u_char &precisionT,
                      const bool &bbox) {

            bytes_t twkb = createHeader(0x05, precisionXY, bbox, false, false, true, false, precisionZ, precisionT);

            if (bbox) {

                int32_t minx = shrink(lines.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(lines.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(lines.front().front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;
                int32_t mint = shrink(lines.front().front().t, static_cast<signed char>(precisionT));
                int32_t maxt = mint;


                for (const auto &line : lines) {

                    for (const auto &location : line) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));
                        int32_t tmpT = shrink(location.t, static_cast<signed char>(precisionT));

                        if (tmpX < minx) minx = tmpX;
                        if (tmpY < miny) miny = tmpY;
                        if (tmpZ < minz) minz = tmpZ;
                        if (tmpT < mint) mint = tmpT;
                        if (tmpX > maxx) maxx = tmpX;
                        if (tmpY > maxy) maxy = tmpY;
                        if (tmpZ > maxz) maxz = tmpZ;
                        if (tmpT > maxt) maxt = tmpT;
                    }
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t bboxT = encodeVarint(encodeZigZag(mint));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));
                bytes_t deltaT = encodeVarint(encodeZigZag(maxt - mint));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ, bboxT, deltaT);

            }

            // Set number of containing points
            bytes_t nlines = encodeVarint(static_cast<uint32_t >(lines.size()));
            append(twkb, nlines);


            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;
            int32_t lastTFull = 0;

            for (const auto &points : lines) {

                bytes_t npoints = encodeVarint(static_cast<uint32_t >(points.size()));
                append(twkb, npoints);

                for (const auto &point : points) {

                    int32_t deltaX = shrink(point.x, precisionXY) - lastXFull;
                    int32_t deltaY = shrink(point.y, precisionXY) - lastYFull;
                    int32_t deltaZ = shrink(point.z, static_cast<signed char>(precisionZ)) - lastZFull;
                    int32_t deltaT = shrink(point.t, static_cast<signed char>(precisionT)) - lastTFull;

                    bytes_t x = encodeVarint(encodeZigZag(deltaX));
                    bytes_t y = encodeVarint(encodeZigZag(deltaY));
                    bytes_t z = encodeVarint(encodeZigZag(deltaZ));
                    bytes_t t = encodeVarint(encodeZigZag(deltaT));

                    append(twkb, x, y, z, t);

                    lastXFull = shrink(point.x, precisionXY);
                    lastYFull = shrink(point.y, precisionXY);
                    lastZFull = shrink(point.z, static_cast<signed char>(precisionZ));
                    lastTFull = shrink(point.t, static_cast<signed char>(precisionT));

                }
            }

            return twkb;
        }

        bytes_t makeMultiLine(const vector<vector<LocationXYZ>> &lines, const signed char &precisionXY, const u_char &precisionZ, const bool &bbox) {

            bytes_t twkb = createHeader(0x05, precisionXY, bbox, false, false, true, false, precisionZ);

            if (bbox) {

                int32_t minx = shrink(lines.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(lines.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(lines.front().front().z, static_cast<signed char>(precisionZ));
                int32_t maxz = minz;


                for (const auto &line : lines) {

                    for (const auto &location : line) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, static_cast<signed char>(precisionZ));

                        if (tmpX < minx) minx = tmpX;
                        if (tmpY < miny) miny = tmpY;
                        if (tmpZ < minz) minz = tmpZ;
                        if (tmpX > maxx) maxx = tmpX;
                        if (tmpY > maxy) maxy = tmpY;
                        if (tmpZ > maxz) maxz = tmpZ;
                    }
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t bboxZ = encodeVarint(encodeZigZag(minz));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));
                bytes_t deltaZ = encodeVarint(encodeZigZag(maxz - minz));

                append(twkb, bboxX, deltaX, bboxY, deltaY, bboxZ, deltaZ);

            }

            // Set number of containing points
            bytes_t nlines = encodeVarint(static_cast<uint32_t >(lines.size()));
            append(twkb, nlines);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;

            for (const auto &points : lines) {

                bytes_t npoints = encodeVarint(static_cast<uint32_t>(points.size()));
                append(twkb, npoints);

                for (const auto &point : points) {

                    int32_t deltaX = shrink(point.x, precisionXY) - lastXFull;
                    int32_t deltaY = shrink(point.y, precisionXY) - lastYFull;
                    int32_t deltaZ = shrink(point.z, static_cast<signed char>(precisionZ)) - lastZFull;

                    bytes_t x = encodeVarint(encodeZigZag(deltaX));
                    bytes_t y = encodeVarint(encodeZigZag(deltaY));
                    bytes_t z = encodeVarint(encodeZigZag(deltaZ));

                    append(twkb, x, y, z);

                    lastXFull = shrink(point.x, precisionXY);
                    lastYFull = shrink(point.y, precisionXY);
                    lastZFull = shrink(point.z, static_cast<signed char>(precisionZ));

                }
            }

            return twkb;
        }

        bytes_t makeMultiLine(const vector<vector<LocationXY>> &lines, const signed char &precisionXY, const bool &bbox) {

            bytes_t twkb = createHeader(0x05, precisionXY, bbox, false, false, false, false);

            if (bbox) {

                int32_t minx = shrink(lines.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(lines.front().front().y, precisionXY);
                int32_t maxy = miny;

                for (const auto &line : lines) {

                    for (const auto &location : line) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);

                        if (tmpX < minx) minx = tmpX;
                        if (tmpY < miny) miny = tmpY;
                        if (tmpX > maxx) maxx = tmpX;
                        if (tmpY > maxy) maxy = tmpY;
                    }
                }

                bytes_t bboxX = encodeVarint(encodeZigZag(minx));
                bytes_t bboxY = encodeVarint(encodeZigZag(miny));
                bytes_t deltaX = encodeVarint(encodeZigZag(maxx - minx));
                bytes_t deltaY = encodeVarint(encodeZigZag(maxy - miny));

                append(twkb, bboxX, deltaX, bboxY, deltaY);

            }

            // Set number of containing points
            bytes_t nlines = encodeVarint(static_cast<uint32_t >(lines.size()));
            append(twkb, nlines);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;

            for (const auto &points : lines) {

                bytes_t npoints = encodeVarint(static_cast<uint32_t >(points.size()));
                append(twkb, npoints);

                for (const auto &point: points) {

                    int32_t deltaX = shrink(point.x, precisionXY) - lastXFull;
                    int32_t deltaY = shrink(point.y, precisionXY) - lastYFull;

                    bytes_t x = encodeVarint(encodeZigZag(deltaX));
                    bytes_t y = encodeVarint(encodeZigZag(deltaY));

                    append(twkb, x, y);

                    lastXFull = shrink(point.x, precisionXY);
                    lastYFull = shrink(point.y, precisionXY);
                }
            }

            return twkb;
        }
    };

};


#endif //TWKB_HEADER_ONLY_TWKB_H
