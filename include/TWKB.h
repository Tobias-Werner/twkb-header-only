//
// Created by Tobias Werner
//

#ifndef TWKB_HEADER_ONLY_TWKB_H
#define TWKB_HEADER_ONLY_TWKB_H

#include <vector>

namespace TWKB {

    using namespace std;

    using bytes_t = vector<unsigned char>;

    struct PosXY {
        double x, y;

        PosXY(double x, double y) : x(x), y(y) {}
    };

    struct PosXYZ {
        double x, y, z;

        PosXYZ(double x, double y, double z) : x(x), y(y), z(z) {}
    };

    struct PosXYZT {
        double x, y, z, t;

        PosXYZT(double x, double y, double z, double t) : x(x), y(y), z(z), t(t) {}
    };


    class GeomFactory {

    private:

        void setTypeAndPrecision(bytes_t &twkb, const unsigned char &type, const signed char &precisionXY) {
            auto &header = twkb[0];
            header |= type;
            header |= encodeZigZag(precisionXY) << 4;
        }

        void setMetadataBits(bytes_t &twkb, const bool &bbox, const bool &sizes, const bool &idlist, const bool &extendedDimensions, const bool &emptyGeom) {
            auto &metadataHeader = twkb[1];
            metadataHeader |= (bbox) ? 0x01 : 0x00;
            metadataHeader |= (sizes) ? 0x02 : 0x00;
            metadataHeader |= (idlist) ? 0x04 : 0x00;
            metadataHeader |= (extendedDimensions) ? 0x08 : 0x00;
            metadataHeader |= (emptyGeom) ? 0x10 : 0x00;
        }

        void addZDimensions(bytes_t &twkb, const unsigned char &precisionZ) {
            unsigned char extendedInformation = 0x01;
            extendedInformation |= (precisionZ << 2);
            twkb.push_back(extendedInformation);
        }

        void addZTDimensions(bytes_t &twkb, const unsigned char &precisionZ, const unsigned char &precisionT) {
            unsigned char extendedInformation = 0x03;
            extendedInformation |= (precisionZ << 2);
            extendedInformation |= (precisionT << 5);
            twkb.push_back(extendedInformation);
        }

        int32_t shrink(const double &value, const signed char &precision) {
            return (precision >= 0) ? round(value * pow10(precision)) : round(value / pow10(-precision));
        }

        bytes_t encode(const double &value, const signed char &precision) {
            return encodeVarint(encodeZigZag(shrink(value, precision)));
        }

        void append(bytes_t &twkb) {}

        template<typename T, typename... Types>
        void append(bytes_t &twkb, const T &var1, const Types... var2) {
            twkb.insert(twkb.end(), var1.begin(), var1.end());
            append(twkb, var2...);
        }

    public:

        static int32_t pow10(unsigned char exponent) {
            static uint32_t pow10[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
            return pow10[exponent];
        }

        static int32_t decodeVarint(const bytes_t &bytes) {

            int32_t result = 0x00;

            for (size_t i = bytes.size(); i > 0; i--) {
                result |= (bytes[i - 1] & 0x7F);
                result <<= (i - 1 > 0) ? 7 : 0;
            }

            return result;
        }

        static bytes_t encodeVarint(int32_t value) {
            bytes_t result;

            do {
                result.push_back(0x80 | (value & 0x7F));
                value >>= 7;
            } while (value > 0x00);

            result.back() &= 0x7F;

            return result;
        }

        static int32_t encodeZigZag(const int32_t &value) {
            return (value << 1) ^ (value >> 31);
        }

        static int32_t decodeZigZag(const uint32_t &value) {
            return -(value & 1) ^ (value >> 1);
        }

        bytes_t makePoint(const PosXY &location, const signed char &precisionXY) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x01, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, false, false);

            // Locations
            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);

            append(twkb, x, y);

            return twkb;
        }

        bytes_t makePoint(const PosXYZ &location, const signed char &precisionXY, const signed char &precisionZ) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x01, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, true, false);

            addZDimensions(twkb, precisionZ);

            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);
            auto z = encode(location.z, precisionZ);

            append(twkb, x, y, z);

            return twkb;
        }

        bytes_t makePoint(const PosXYZT &location, const signed char &precisionXY, const signed char &precisionZ, const signed char &precisionT) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x01, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, false, false, false, true, false);

            addZTDimensions(twkb, precisionZ, precisionT);

            auto x = encode(location.x, precisionXY);
            auto y = encode(location.y, precisionXY);
            auto z = encode(location.z, precisionZ);
            auto t = encode(location.t, precisionT);

            append(twkb, x, y, z, t);

            return twkb;
        }

        bytes_t makeLine(const vector<PosXY> &locations, const signed char &precisionXY, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x02, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, false, false);

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
            bytes_t npoints = encodeVarint(locations.size());

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

        bytes_t makeLine(const vector<PosXYZ> &locations, const signed char &precisionXY, const signed char &precisionZ, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x02, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            // Add extendet dimension information
            addZDimensions(twkb, precisionZ);

            // Add bounding box - if selected
            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, precisionZ);
                int32_t maxz = minz;

                for (const auto &location : locations) {
                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, precisionZ);

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
            bytes_t npoints = encodeVarint(locations.size());
            twkb.insert(twkb.end(), npoints.begin(), npoints.end());

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, precisionZ);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));

            append(twkb, x, y, z);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, precisionZ) - zShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));

                append(twkb, x, y, z);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, precisionZ);

            }

            return twkb;
        }

        bytes_t makeLine(const vector<PosXYZT> &locations, const signed char &precisionXY, const signed char &precisionZ, const signed char &precisionT, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x02, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            // Add extendet dimension information
            addZTDimensions(twkb, precisionZ, precisionT);

            // Add bounding box - if selected
            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, precisionZ);
                int32_t maxz = minz;
                int32_t mint = shrink(locations.front().t, precisionT);
                int32_t maxt = mint;

                for (const auto &location : locations) {
                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, precisionZ);
                    int32_t tmpT = shrink(location.t, precisionT);

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
            bytes_t npoints = encodeVarint(locations.size());
            twkb.insert(twkb.end(), npoints.begin(), npoints.end());

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, precisionZ);
            int32_t tShrinked = shrink(locations.front().t, precisionT);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));
            bytes_t t = encodeVarint(encodeZigZag(tShrinked));

            append(twkb, x, y, z, t);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, precisionZ) - zShrinked;
                int32_t deltaT = shrink(locations[i].t, precisionT) - tShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));
                t = encodeVarint(encodeZigZag(deltaT));

                append(twkb, x, y, z, t);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, precisionZ);
                tShrinked = shrink(locations[i].t, precisionT);

            }

            return twkb;
        }

        bytes_t makePolygon(const vector<vector<PosXY>> &rings, const signed char &precisionXY, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x03, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, false, false);


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
            bytes_t nrings = encodeVarint(rings.size());
            append(twkb, nrings);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;

            for (size_t i = 0; i < rings.size(); i++) {
                auto &pointsInRing = rings[i];

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(pointsInRing.size());
                append(twkb, npoints);

                for (size_t j = 0; j < pointsInRing.size(); j++) {

                    int32_t currentXFull = shrink(pointsInRing[j].x, precisionXY);
                    int32_t currentYFull = shrink(pointsInRing[j].y, precisionXY);

                    int32_t deltaX = currentXFull - lastXFull;
                    int32_t deltaY = currentYFull - lastYFull;

                    lastXFull = currentXFull;
                    lastYFull = currentYFull;

                    bytes_t deltaXAsZigZag = encodeVarint(encodeZigZag(deltaX));
                    bytes_t deltaYAsZigZag = encodeVarint(encodeZigZag(deltaY));

                    append(twkb, deltaXAsZigZag, deltaYAsZigZag);

                }
            }


            return twkb;
        }

        bytes_t makePolygon(const vector<vector<PosXYZ>> &rings, const signed char &precisionXY, const signed char &precisionZ, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x03, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            addZDimensions(twkb, precisionZ);

            if (bbox) {

                int32_t minx = shrink(rings.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(rings.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(rings.front().front().z, precisionZ);
                int32_t maxz = minz;

                for (const auto &ring : rings) {

                    for (const auto &location : ring) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, precisionZ);

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
            bytes_t nrings = encodeVarint(rings.size());
            append(twkb, nrings);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;

            for (size_t i = 0; i < rings.size(); i++) {
                const auto &pointsInRing = rings[i];

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(pointsInRing.size());
                append(twkb, npoints);

                for (size_t j = 0; j < pointsInRing.size(); j++) {

                    int32_t currentXFull = shrink(pointsInRing[j].x, precisionXY);
                    int32_t currentYFull = shrink(pointsInRing[j].y, precisionXY);
                    int32_t currentZFull = shrink(pointsInRing[j].z, precisionZ);

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

        bytes_t makePolygon(const vector<vector<PosXYZT>> &rings, const signed char &precisionXY, const signed char &precisionZ, const signed char &precisionT, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x03, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            addZTDimensions(twkb, precisionZ, precisionT);

            if (bbox) {

                int32_t minx = shrink(rings.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(rings.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(rings.front().front().z, precisionZ);
                int32_t maxz = minz;
                int32_t mint = shrink(rings.front().front().t, precisionT);
                int32_t maxt = mint;

                for (const auto &ring : rings) {

                    for (const auto &location : ring) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, precisionZ);
                        int32_t tmpT = shrink(location.t, precisionT);

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
            bytes_t nrings = encodeVarint(rings.size());
            append(twkb, nrings);


            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;
            int32_t lastTFull = 0;

            for (size_t i = 0; i < rings.size(); i++) {
                auto &pointsInRing = rings[i];

                // Set number of containing pointsInRing
                bytes_t npoints = encodeVarint(pointsInRing.size());
                append(twkb, npoints);

                for (size_t j = 0; j < pointsInRing.size(); j++) {

                    int32_t currentXFull = shrink(pointsInRing[j].x, precisionXY);
                    int32_t currentYFull = shrink(pointsInRing[j].y, precisionXY);
                    int32_t currentZFull = shrink(pointsInRing[j].z, precisionZ);
                    int32_t currentTFull = shrink(pointsInRing[j].t, precisionT);

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

        bytes_t makeMultiPoint(const vector<PosXY> &locations, const signed char &precisionXY, const bool &bbox) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x04, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, false, false);

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
            bytes_t npoints = encodeVarint(locations.size());
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

        bytes_t makeMultiPoint(const vector<PosXYZ> &locations, const signed char &precisionXY, const signed char &precisionZ, const bool &bbox) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x04, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            addZDimensions(twkb, precisionZ);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, precisionZ);
                int32_t maxz = minz;

                for (const auto &location : locations) {

                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, precisionZ);

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
            bytes_t npoints = encodeVarint(locations.size());
            append(twkb, npoints);

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, precisionZ);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));

            append(twkb, x, y, z);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, precisionZ) - zShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));

                append(twkb, x, y, z);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, precisionZ);

            }


            return twkb;
        }

        bytes_t makeMultiPoint(const vector<PosXYZT> &locations, const signed char &precisionXY, const signed char &precisionZ, const signed char &precisionT, const bool &bbox) {

            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x04, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            addZTDimensions(twkb, precisionZ, precisionT);

            if (bbox) {

                int32_t minx = shrink(locations.front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(locations.front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(locations.front().z, precisionZ);
                int32_t maxz = minz;
                int32_t mint = shrink(locations.front().t, precisionT);
                int32_t maxt = mint;

                for (const auto &location : locations) {

                    int32_t tmpX = shrink(location.x, precisionXY);
                    int32_t tmpY = shrink(location.y, precisionXY);
                    int32_t tmpZ = shrink(location.z, precisionZ);
                    int32_t tmpT = shrink(location.t, precisionT);

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
            bytes_t npoints = encodeVarint(locations.size());
            append(twkb, npoints);

            int32_t xShrinked = shrink(locations.front().x, precisionXY);
            int32_t yShrinked = shrink(locations.front().y, precisionXY);
            int32_t zShrinked = shrink(locations.front().z, precisionZ);
            int32_t tShrinked = shrink(locations.front().t, precisionT);

            bytes_t x = encodeVarint(encodeZigZag(xShrinked));
            bytes_t y = encodeVarint(encodeZigZag(yShrinked));
            bytes_t z = encodeVarint(encodeZigZag(zShrinked));
            bytes_t t = encodeVarint(encodeZigZag(tShrinked));

            append(twkb, x, y, z, t);

            for (size_t i = 1; i < locations.size(); i++) {

                int32_t deltaX = shrink(locations[i].x, precisionXY) - xShrinked;
                int32_t deltaY = shrink(locations[i].y, precisionXY) - yShrinked;
                int32_t deltaZ = shrink(locations[i].z, precisionZ) - zShrinked;
                int32_t deltaT = shrink(locations[i].t, precisionT) - tShrinked;

                x = encodeVarint(encodeZigZag(deltaX));
                y = encodeVarint(encodeZigZag(deltaY));
                z = encodeVarint(encodeZigZag(deltaZ));
                t = encodeVarint(encodeZigZag(deltaT));

                append(twkb, x, y, z, t);

                xShrinked = shrink(locations[i].x, precisionXY);
                yShrinked = shrink(locations[i].y, precisionXY);
                zShrinked = shrink(locations[i].z, precisionZ);
                tShrinked = shrink(locations[i].t, precisionT);

            }


            return twkb;
        }

        bytes_t
        makeMultiLine(const vector<vector<PosXYZT>> &lines, const signed char &precisionXY, const signed char &precisionZ, const signed char &precisionT, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x05, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            // Add extendet dimension information
            addZTDimensions(twkb, precisionZ, precisionT);

            // Add bounding box - if selected
            if (bbox) {

                int32_t minx = shrink(lines.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(lines.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(lines.front().front().z, precisionZ);
                int32_t maxz = minz;
                int32_t mint = shrink(lines.front().front().t, precisionT);
                int32_t maxt = mint;


                for (const auto &line : lines) {

                    for (const auto &location : line) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, precisionZ);
                        int32_t tmpT = shrink(location.t, precisionT);

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
            bytes_t nlines = encodeVarint(lines.size());
            append(twkb, nlines);


            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;
            int32_t lastTFull = 0;

            for (size_t i = 0; i < lines.size(); i++) {

                const auto &points = lines[i];

                bytes_t npoints = encodeVarint(points.size());
                append(twkb, npoints);

                for (size_t j = 0; j < points.size(); j++) {

                    const auto &point = points[j];

                    int32_t deltaX = shrink(point.x, precisionXY) - lastXFull;
                    int32_t deltaY = shrink(point.y, precisionXY) - lastYFull;
                    int32_t deltaZ = shrink(point.z, precisionZ) - lastZFull;
                    int32_t deltaT = shrink(point.t, precisionT) - lastTFull;

                    bytes_t x = encodeVarint(encodeZigZag(deltaX));
                    bytes_t y = encodeVarint(encodeZigZag(deltaY));
                    bytes_t z = encodeVarint(encodeZigZag(deltaZ));
                    bytes_t t = encodeVarint(encodeZigZag(deltaT));

                    append(twkb, x, y, z, t);

                    lastXFull = shrink(point.x, precisionXY);
                    lastYFull = shrink(point.y, precisionXY);
                    lastZFull = shrink(point.z, precisionZ);
                    lastTFull = shrink(point.t, precisionT);

                }
            }

            return twkb;
        }

        bytes_t makeMultiLine(const vector<vector<PosXYZ>> &lines, const signed char &precisionXY, const signed char &precisionZ, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x05, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, true, false);

            // Add extendet dimension information
            addZDimensions(twkb, precisionZ);

            // Add bounding box - if selected
            if (bbox) {

                int32_t minx = shrink(lines.front().front().x, precisionXY);
                int32_t maxx = minx;
                int32_t miny = shrink(lines.front().front().y, precisionXY);
                int32_t maxy = miny;
                int32_t minz = shrink(lines.front().front().z, precisionZ);
                int32_t maxz = minz;


                for (const auto &line : lines) {

                    for (const auto &location : line) {
                        int32_t tmpX = shrink(location.x, precisionXY);
                        int32_t tmpY = shrink(location.y, precisionXY);
                        int32_t tmpZ = shrink(location.z, precisionZ);

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
            bytes_t nlines = encodeVarint(lines.size());
            append(twkb, nlines);


            int32_t lastXFull = 0;
            int32_t lastYFull = 0;
            int32_t lastZFull = 0;

            for (size_t i = 0; i < lines.size(); i++) {

                const auto &points = lines[i];

                bytes_t npoints = encodeVarint(points.size());
                append(twkb, npoints);

                for (size_t j = 0; j < points.size(); j++) {

                    const auto &point = points[j];

                    int32_t deltaX = shrink(point.x, precisionXY) - lastXFull;
                    int32_t deltaY = shrink(point.y, precisionXY) - lastYFull;
                    int32_t deltaZ = shrink(point.z, precisionZ) - lastZFull;

                    bytes_t x = encodeVarint(encodeZigZag(deltaX));
                    bytes_t y = encodeVarint(encodeZigZag(deltaY));
                    bytes_t z = encodeVarint(encodeZigZag(deltaZ));

                    append(twkb, x, y, z);

                    lastXFull = shrink(point.x, precisionXY);
                    lastYFull = shrink(point.y, precisionXY);
                    lastZFull = shrink(point.z, precisionZ);

                }
            }

            return twkb;
        }

        bytes_t makeMultiLine(const vector<vector<PosXY>> &lines, const signed char &precisionXY, const bool &bbox) {
            bytes_t twkb({0x00, 0x00});

            // Set type and precision
            setTypeAndPrecision(twkb, 0x05, precisionXY);

            // Prepare metadata header
            // bbox | size | idlist | extended dimensions | empty geom
            setMetadataBits(twkb, bbox, false, false, false, false);

            // Add bounding box - if selected
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
            bytes_t nlines = encodeVarint(lines.size());
            append(twkb, nlines);

            int32_t lastXFull = 0;
            int32_t lastYFull = 0;

            for (size_t i = 0; i < lines.size(); i++) {

                const auto &points = lines[i];

                bytes_t npoints = encodeVarint(points.size());
                append(twkb, npoints);

                for (size_t j = 0; j < points.size(); j++) {

                    const auto &point = points[j];

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
