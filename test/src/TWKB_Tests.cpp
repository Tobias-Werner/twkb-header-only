//
// Created by Tobias Werner
//
// TWKB target bytes are generated by PostgreSQL/Postgis

#define CATCH_CONFIG_MAIN

#include <stdio.h>
#include <iostream>
#include <catch.hpp>
#include "TWKB.h"

using namespace TWKB;
using namespace std;


TEST_CASE("Binary encoding") {

    SECTION("Varint encoding") {

        vector<unsigned int> testValues({0, 300, 12434, 324345345});

        for (auto &testValue : testValues) {
            bytes_t bytes = GeomFactory::encodeVarint(testValue);
            CHECK(GeomFactory::decodeVarint(bytes) == testValue);
        }

    }

    SECTION("ZigZag encoding") {
        vector<int> testValues({0, 300, -300, 2545645, -2545645});

        for (auto &testValue : testValues) {
            unsigned int encoded = GeomFactory::encodeZigZag(testValue);
            CHECK(GeomFactory::decodeZigZag(encoded) == testValue);
        }
    }

}


TEST_CASE("Creating point geometries") {

    SECTION("Creating point XY") {

        GeomFactory factory;

        PosXY pos(7.625752, 53.942254);
        bytes_t twkb = factory.makePoint(pos, 6);

        // SELECT ST_AsTWKB('POINT(7.625752 53.942254)'::geometry, 6) as binary;
        bytes_t targetTwkb = bytes_t({0xC1, 0x00, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33});

        CHECK(twkb == targetTwkb);
    }


    SECTION("Creating point XYZ") {

        GeomFactory factory;

        PosXYZ pos(7.625752, 53.942254, 10.175);
        bytes_t twkb = factory.makePoint(pos, 6, 3);

        // SELECT ST_AsTWKB('POINT(7.625752 53.942254 10.175)'::geometry, 6, 3) as binary;
        bytes_t targetTwkb = bytes_t({0xC1, 0x08, 0x0D, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01});

        CHECK(twkb == targetTwkb);
    }

    SECTION("Creating point XYZT") {

        GeomFactory factory;

        PosXYZT pos(7.625752, 53.942254, 10.175, 164.0);
        bytes_t twkb = factory.makePoint(pos, 6, 3, 0);

        // SELECT ST_AsTWKB('POINT(7.625752 53.942254 10.175 164.0)'::geometry, 6, 3, 0) as binary;
        bytes_t targetTwkb = bytes_t({0xC1, 0x08, 0x0F, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0xC8, 0x02});

        CHECK(twkb == targetTwkb);
    }
}

TEST_CASE("Creating line geometries") {

    SECTION("Creating line XY") {
        GeomFactory factory;

        PosXY pos1(7.625752, 53.942254);
        PosXY pos2(7.615752, 53.932254);
        PosXY pos3(7.532752, 53.915354);

        vector<PosXY> locations({pos1, pos2, pos3});

        bytes_t twkb = factory.makeLine(locations, 6, false);

        // SELECT ST_AsTWKB('LINESTRING(7.625752 53.942254, 7.615752 53.932254, 7.532752 53.915354)'::geometry, 6) as binary;
        bytes_t targetTwkb = bytes_t({0xC2, 0x00, 0x03, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating line XY - with bbox") {
        GeomFactory factory;

        PosXY pos1(7.625752, 53.942254);
        PosXY pos2(7.615752, 53.932254);
        PosXY pos3(7.532752, 53.915354);

        vector<PosXY> locations({pos1, pos2, pos3});

        bytes_t twkb = factory.makeLine(locations, 6, true);

        // SELECT ST_AsTWKB('LINESTRING(7.625752 53.942254, 7.615752 53.932254, 7.532752 53.915354)'::geometry, 6) as binary;
        bytes_t targetTwkb = bytes_t(
                {0xC2, 0x01, 0xA0, 0xC3, 0x97, 0x07, 0x90, 0xAD, 0x0B, 0xB4, 0xBB, 0xB5, 0x33, 0xA8, 0xA4, 0x03, 0x03, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0x9F, 0x9C,
                 0x01, 0x9F, 0x9C, 0x01, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating line XYZ") {
        GeomFactory factory;

        PosXYZ pos1(7.625752, 53.942254, 10.175);
        PosXYZ pos2(7.615752, 53.932254, 10.231);
        PosXYZ pos3(7.532752, 53.915354, 10.335);

        vector<PosXYZ> locations({pos1, pos2, pos3});

        bytes_t twkb = factory.makeLine(locations, 6, 3, false);

        // SELECT ST_AsTWKB('LINESTRING(7.625752 53.942254 10.175, 7.615752 53.932254 10.231, 7.532752 53.915354 10.335)'::geometry, 6, 3) as binary;
        bytes_t targetTwkb = bytes_t(
                {0xC2, 0x08, 0x0D, 0x03, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0xEF, 0x90, 0x0A, 0x87, 0x88,
                 0x02, 0xD0, 0x01});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating line XYZ - with bbox") {
        GeomFactory factory;

        PosXYZ pos1(7.625752, 53.942254, 10.175);
        PosXYZ pos2(7.615752, 53.932254, 10.231);
        PosXYZ pos3(7.532752, 53.915354, 10.335);

        vector<PosXYZ> locations({pos1, pos2, pos3});

        bytes_t twkb = factory.makeLine(locations, 6, 3, true);

        // SELECT ST_AsTWKB('LINESTRING(7.625752 53.942254 10.175, 7.615752 53.932254 10.231, 7.532752 53.915354 10.335)'::geometry, 6, 3, 0, false, true) as binary;
        bytes_t targetTwkb = bytes_t(
                {0xC2, 0x09, 0x0D, 0xA0, 0xC3, 0x97, 0x07, 0x90, 0xAD, 0x0B, 0xB4, 0xBB, 0xB5, 0x33, 0xA8, 0xA4, 0x03, 0xFE, 0x9E, 0x01, 0xC0, 0x02, 0x03, 0xB0, 0xF0, 0xA2, 0x07,
                 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02, 0xD0, 0x01});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating line XYZT") {
        GeomFactory factory;

        PosXYZT pos1(7.625752, 53.942254, 10.175, 164.0);
        PosXYZT pos2(7.615752, 53.932254, 10.231, 165.0);
        PosXYZT pos3(7.532752, 53.915354, 10.335, 166.0);

        vector<PosXYZT> locations({pos1, pos2, pos3});

        bytes_t twkb = factory.makeLine(locations, 6, 3, 0, false);

        // SELECT ST_AsTWKB('LINESTRING(7.625752 53.942254 10.175 164.0, 7.615752 53.932254 10.231 165.0, 7.532752 53.915354 10.335 166.0)'::geometry, 6, 3, 0) as binary;
        bytes_t targetTwkb = bytes_t(
                {0xC2, 0x08, 0x0F, 0x03, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0xC8, 0x02, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0x02, 0xEF, 0x90,
                 0x0A, 0x87, 0x88, 0x02, 0xD0, 0x01, 0x02});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating line XYZT - with bbox") {
        GeomFactory factory;

        PosXYZT pos1(7.625752, 53.942254, 10.175, 164.0);
        PosXYZT pos2(7.615752, 53.932254, 10.231, 165.0);
        PosXYZT pos3(7.532752, 53.915354, 10.335, 166.0);

        vector<PosXYZT> locations({pos1, pos2, pos3});

        bytes_t twkb = factory.makeLine(locations, 6, 3, 0, true);

        // SELECT ST_AsTWKB('LINESTRING(7.625752 53.942254 10.175 164.0, 7.615752 53.932254 10.231 165.0, 7.532752 53.915354 10.335 166.0)'::geometry, 6, 3, 0, false, true) as binary;
        bytes_t targetTwkb = bytes_t(
                {0xC2, 0x09, 0x0F, 0xA0, 0xC3, 0x97, 0x07, 0x90, 0xAD, 0x0B, 0xB4, 0xBB, 0xB5, 0x33, 0xA8, 0xA4, 0x03, 0xFE, 0x9E, 0x01, 0xC0, 0x02, 0xC8, 0x02, 0x04, 0x03, 0xB0,
                 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0xC8, 0x02, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0x02, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02, 0xD0,
                 0x01, 0x02});

        CHECK(twkb == targetTwkb);

    }
}

TEST_CASE("Creating polygon geometries") {

    SECTION("Creating polygon XY") {
        GeomFactory factory;

        PosXY pos1(7.625752, 53.942254);
        PosXY pos2(7.615752, 53.932254);
        PosXY pos3(7.532752, 53.915354);

        vector<PosXY> exteriorRing({pos1, pos2, pos3, pos1});
        vector<vector<PosXY>> polygon({exteriorRing});

        PosXY pos4(7.60679, 53.93277);
        PosXY pos5(7.61170, 53.93686);
        PosXY pos6(7.61577, 53.93547);

        vector<PosXY> interiorRing({pos4, pos5, pos6, pos4});
        polygon.push_back(interiorRing);

        bytes_t twkb = factory.makePolygon(polygon, 6, false);

        // NOTE: TWKB doesn't store duplicated coordinates for closing rings. Instead of this, they are implicity closed.
        // SELECT ST_AsTWKB('POLYGON((7.625752 53.942254, 7.615752 53.932254, 7.532752 53.915354, 7.625752 53.942254), (7.60679 53.93277, 7.61170 53.93686, 7.61577 53.93547, 7.60679 53.93277))'::geometry, 6) as binary;

        bytes_t targetTwkb = bytes_t(
                {0xC3, 0x00, 0x02, 0x04, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02, 0x90, 0xAD, 0x0B,
                 0xA8, 0xA4, 0x03, 0x04, 0xA3, 0xA8, 0x02, 0x97, 0x94, 0x01, 0xDC, 0x4C, 0xF4, 0x3F, 0xCC, 0x3F, 0xDB, 0x15, 0xA7, 0x8C, 0x01, 0x97, 0x2A});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating polygon XY - with bbox") {
        GeomFactory factory;

        PosXY pos1(7.625752, 53.942254);
        PosXY pos2(7.615752, 53.932254);
        PosXY pos3(7.532752, 53.915354);

        vector<PosXY> exteriorRing({pos1, pos2, pos3, pos1});
        vector<vector<PosXY>> polygon({exteriorRing});

        PosXY pos4(7.60679, 53.93277);
        PosXY pos5(7.61170, 53.93686);
        PosXY pos6(7.61577, 53.93547);

        vector<PosXY> interiorRing({pos4, pos5, pos6, pos4});
        polygon.push_back(interiorRing);

        bytes_t twkb = factory.makePolygon(polygon, 6, true);

        // NOTE: TWKB doesn't store duplicated coordinates for closing rings. Instead of this, they are implicity closed.
        // SELECT ST_AsTWKB('POLYGON((7.625752 53.942254, 7.615752 53.932254, 7.532752 53.915354, 7.625752 53.942254), (7.60679 53.93277, 7.61170 53.93686, 7.61577 53.93547, 7.60679 53.93277))'::geometry, 6, 0, 0, false, true) as binary;

        bytes_t targetTwkb = bytes_t(
                {0xC3, 0x01, 0xA0, 0xC3, 0x97, 0x07, 0x90, 0xAD, 0x0B, 0xB4, 0xBB, 0xB5, 0x33, 0xA8, 0xA4, 0x03, 0x02, 0x04, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0x9F,
                 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02, 0x90, 0xAD, 0x0B, 0xA8, 0xA4, 0x03, 0x04, 0xA3, 0xA8, 0x02, 0x97, 0x94, 0x01, 0xDC, 0x4C, 0xF4,
                 0x3F, 0xCC, 0x3F, 0xDB, 0x15, 0xA7, 0x8C, 0x01, 0x97, 0x2A});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating polygon XYZ") {
        GeomFactory factory;


        PosXYZ pos1(7.625752, 53.942254, 10.175);
        PosXYZ pos2(7.615752, 53.932254, 10.231);
        PosXYZ pos3(7.532752, 53.915354, 10.335);

        vector<PosXYZ> exteriorRing({pos1, pos2, pos3, pos1});
        vector<vector<PosXYZ>> polygon({exteriorRing});

        PosXYZ pos4(7.60679, 53.93277, 10.155);
        PosXYZ pos5(7.61170, 53.93686, 10.158);
        PosXYZ pos6(7.61577, 53.93547, 10.651);

        vector<PosXYZ> interiorRing({pos4, pos5, pos6, pos4});
        polygon.push_back(interiorRing);

        bytes_t twkb = factory.makePolygon(polygon, 6, 3, false);

        // NOTE: TWKB doesn't store duplicated coordinates for closing rings. Instead of this, they are implicity closed.
        // SELECT ST_AsTWKB('POLYGON((7.625752 53.942254 10.175, 7.615752 53.932254 10.231, 7.532752 53.915354 10.335, 7.625752 53.942254 10.175), (7.60679 53.93277 10.155, 7.61170 53.93686 10.158, 7.61577 53.93547 10.651, 7.60679 53.93277 10.155))'::geometry, 6, 3, 0, false, false) as binary;

        bytes_t targetTwkb = bytes_t(
                {0xC3, 0x08, 0x0D, 0x02, 0x04, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0xEF, 0x90, 0x0A, 0x87,
                 0x88, 0x02, 0xD0, 0x01, 0x90, 0xAD, 0x0B, 0xA8, 0xA4, 0x03, 0xBF, 0x02, 0x04, 0xA3, 0xA8, 0x02, 0x97, 0x94, 0x01, 0x27, 0xDC, 0x4C, 0xF4, 0x3F, 0x06, 0xCC, 0x3F,
                 0xDB, 0x15, 0xDA, 0x07, 0xA7, 0x8C, 0x01, 0x97, 0x2A, 0xDF, 0x07});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating polygon XYZ - with bbox") {
        GeomFactory factory;

        PosXYZ pos1(7.625752, 53.942254, 10.175);
        PosXYZ pos2(7.615752, 53.932254, 10.231);
        PosXYZ pos3(7.532752, 53.915354, 10.335);

        vector<PosXYZ> exteriorRing({pos1, pos2, pos3, pos1});
        vector<vector<PosXYZ>> polygon({exteriorRing});

        PosXYZ pos4(7.60679, 53.93277, 10.155);
        PosXYZ pos5(7.61170, 53.93686, 10.158);
        PosXYZ pos6(7.61577, 53.93547, 10.651);

        vector<PosXYZ> interiorRing({pos4, pos5, pos6, pos4});
        polygon.push_back(interiorRing);

        bytes_t twkb = factory.makePolygon(polygon, 6, 3, true);

        // NOTE: TWKB doesn't store duplicated coordinates for closing rings. Instead of this, they are implicity closed.
        // SELECT ST_AsTWKB('POLYGON((7.625752 53.942254 10.175, 7.615752 53.932254 10.231, 7.532752 53.915354 10.335, 7.625752 53.942254 10.175), (7.60679 53.93277 10.155, 7.61170 53.93686 10.158, 7.61577 53.93547 10.651, 7.60679 53.93277 10.155))'::geometry, 6, 3, 0, false, true) as binary;

        bytes_t targetTwkb = bytes_t(
                {0xC3, 0x09, 0x0D, 0xA0, 0xC3, 0x97, 0x07, 0x90, 0xAD, 0x0B, 0xB4, 0xBB, 0xB5, 0x33, 0xA8, 0xA4, 0x03, 0xD6, 0x9E, 0x01, 0xE0, 0x07, 0x02, 0x04, 0xB0, 0xF0, 0xA2,
                 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02, 0xD0, 0x01, 0x90, 0xAD, 0x0B, 0xA8,
                 0xA4, 0x03, 0xBF, 0x02, 0x04, 0xA3, 0xA8, 0x02, 0x97, 0x94, 0x01, 0x27, 0xDC, 0x4C, 0xF4, 0x3F, 0x06, 0xCC, 0x3F, 0xDB, 0x15, 0xDA, 0x07, 0xA7, 0x8C, 0x01, 0x97,
                 0x2A, 0xDF, 0x07});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating polygon XYZT") {
        GeomFactory factory;


        PosXYZT pos1(7.625752, 53.942254, 10.175, 164.0);
        PosXYZT pos2(7.615752, 53.932254, 10.231, 165.0);
        PosXYZT pos3(7.532752, 53.915354, 10.335, 166.0);

        vector<PosXYZT> exteriorRing({pos1, pos2, pos3, pos1});
        vector<vector<PosXYZT>> polygon({exteriorRing});

        PosXYZT pos4(7.60679, 53.93277, 10.155, 171.0);
        PosXYZT pos5(7.61170, 53.93686, 10.158, 172.0);
        PosXYZT pos6(7.61577, 53.93547, 10.651, 173.0);

        vector<PosXYZT> interiorRing({pos4, pos5, pos6, pos4});
        polygon.push_back(interiorRing);

        bytes_t twkb = factory.makePolygon(polygon, 6, 3, 0, false);

        // NOTE: TWKB doesn't store duplicated coordinates for closing rings. Instead of this, they are implicity closed.
        // ST_AsTWKB('POLYGON((7.625752 53.942254 10.175 164.0, 7.615752 53.932254 10.231 165.0, 7.532752 53.915354 10.335 166.0, 7.625752 53.942254 10.175 164.0), (7.60679 53.93277 10.155 171.0, 7.61170 53.93686 10.158 172.0, 7.61577 53.93547 10.651 173.0, 7.60679 53.93277 10.155 171.0))'::geometry, 6, 3, 0, false, false) as binary;

        bytes_t targetTwkb = bytes_t(
                {0xC3, 0x08, 0x0F, 0x02, 0x04, 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0xC8, 0x02, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0x02, 0xEF,
                 0x90, 0x0A, 0x87, 0x88, 0x02, 0xD0, 0x01, 0x02, 0x90, 0xAD, 0x0B, 0xA8, 0xA4, 0x03, 0xBF, 0x02, 0x03, 0x04, 0xA3, 0xA8, 0x02, 0x97, 0x94, 0x01, 0x27, 0x0E, 0xDC,
                 0x4C, 0xF4, 0x3F, 0x06, 0x02, 0xCC, 0x3F, 0xDB, 0x15, 0xDA, 0x07, 0x02, 0xA7, 0x8C, 0x01, 0x97, 0x2A, 0xDF, 0x07, 0x03});

        CHECK(twkb == targetTwkb);

    }

    SECTION("Creating polygon XYZT - with bbox") {
        GeomFactory factory;

        PosXYZT pos1(7.625752, 53.942254, 10.175, 164.0);
        PosXYZT pos2(7.615752, 53.932254, 10.231, 165.0);
        PosXYZT pos3(7.532752, 53.915354, 10.335, 166.0);

        vector<PosXYZT> exteriorRing({pos1, pos2, pos3, pos1});
        vector<vector<PosXYZT>> polygon({exteriorRing});

        PosXYZT pos4(7.60679, 53.93277, 10.155, 171.0);
        PosXYZT pos5(7.61170, 53.93686, 10.158, 172.0);
        PosXYZT pos6(7.61577, 53.93547, 10.651, 173.0);

        vector<PosXYZT> interiorRing({pos4, pos5, pos6, pos4});
        polygon.push_back(interiorRing);

        bytes_t twkb = factory.makePolygon(polygon, 6, 3, 0, true);

        // NOTE: TWKB doesn't store duplicated coordinates for closing rings. Instead of this, they are implicity closed.
        // ST_AsTWKB('POLYGON((7.625752 53.942254 10.175 164.0, 7.615752 53.932254 10.231 165.0, 7.532752 53.915354 10.335 166.0, 7.625752 53.942254 10.175 164.0), (7.60679 53.93277 10.155 171.0, 7.61170 53.93686 10.158 172.0, 7.61577 53.93547 10.651 173.0, 7.60679 53.93277 10.155 171.0))'::geometry, 6, 3, 0, false, true) as binary;

        bytes_t targetTwkb = bytes_t(
                {0xC3, 0x09, 0x0F, 0xA0, 0xC3, 0x97, 0x07, 0x90, 0xAD, 0x0B, 0xB4, 0xBB, 0xB5, 0x33, 0xA8, 0xA4, 0x03, 0xD6, 0x9E, 0x01, 0xE0, 0x07, 0xC8, 0x02, 0x12, 0x02, 0x04,
                 0xB0, 0xF0, 0xA2, 0x07, 0xDC, 0xDF, 0xB8, 0x33, 0xFE, 0x9E, 0x01, 0xC8, 0x02, 0x9F, 0x9C, 0x01, 0x9F, 0x9C, 0x01, 0x70, 0x02, 0xEF, 0x90, 0x0A, 0x87, 0x88, 0x02,
                 0xD0, 0x01, 0x02, 0x90, 0xAD, 0x0B, 0xA8, 0xA4, 0x03, 0xBF, 0x02, 0x03, 0x04, 0xA3, 0xA8, 0x02, 0x97, 0x94, 0x01, 0x27, 0x0E, 0xDC, 0x4C, 0xF4, 0x3F, 0x06, 0x02,
                 0xCC, 0x3F, 0xDB, 0x15, 0xDA, 0x07, 0x02, 0xA7, 0x8C, 0x01, 0x97, 0x2A, 0xDF, 0x07, 0x03});

        CHECK(twkb == targetTwkb);

    }

}