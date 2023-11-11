#include "convert.h"
#include <proj/crs.hpp>
#include <proj/io.hpp>
#include <proj/util.hpp>
#include <vector>
#include <cmath>


using namespace geovi;

bool CoordinateSystemConverter::convert(CoordinateSystemType srcCRS,CoordinateSystemType targetCRS,Point2& point){
    if( srcCRS == CoordinateSystemType::WGS84 && targetCRS == CoordinateSystemType::UTM){
        return convertBetweenWGS84AndUTM(point);
    }
    if( srcCRS == CoordinateSystemType::UTM && targetCRS == CoordinateSystemType::WGS84){
        return convertBetweenUTMAndWGS84(point);
    }
    return false;
}

bool CoordinateSystemConverter::convertBetweenUTMAndWGS84(Point2& point){
    
   // auto dbContext = DatabaseContext::create();
    auto authFactory = AuthorityFactory::create(dbContext,std::string());
    auto coord_op_ctxt = CoordinateOperationContext::create(authFactory,nullptr,0.0);
    auto authFactoryEPSG = AuthorityFactory::create(dbContext,"EPSG");
    auto sourceCRS = authFactoryEPSG->createCoordinateReferenceSystem(std::to_string(int(band)));
    auto targetCRS = authFactoryEPSG->createCoordinateReferenceSystem("4326");
    auto list = CoordinateOperationFactory::create()->createOperations(
        sourceCRS, targetCRS, coord_op_ctxt);
    if(list.empty()){
        return false;
    }
    PJ_CONTEXT *ctx = proj_context_create();
    auto transformer = list[0]->coordinateTransformer(ctx);
    PJ_COORD c = {{
        point.x,
        point.y,
        0.0,
        HUGE_VAL
    }};
    c = transformer->transform(c);
    point.x = c.v[0];
    point.y = c.v[1];
    proj_context_destroy(ctx);

    return true;
}

bool CoordinateSystemConverter::convertBetweenWGS84AndUTM(Point2& point){
    
    //auto dbContext = DatabaseContext::create();
    auto authFactory = AuthorityFactory::create(dbContext,std::string());
    auto coord_op_ctxt = CoordinateOperationContext::create(authFactory,nullptr,0.0);
    auto authFactoryEPSG = AuthorityFactory::create(dbContext,"EPSG");
    auto sourceCRS = authFactoryEPSG->createCoordinateReferenceSystem("4326");
    auto targetCRS = authFactoryEPSG->createCoordinateReferenceSystem(std::to_string(int(band)));
     auto list = CoordinateOperationFactory::create()->createOperations(
        sourceCRS, targetCRS, coord_op_ctxt);
    if(list.empty()){
        return false;
    }
    PJ_CONTEXT *ctx = proj_context_create();
    auto transformer = list[0]->coordinateTransformer(ctx);
    PJ_COORD c = {{
        point.x,
        point.y,
        0.0,
        HUGE_VAL
    }};
    c = transformer->transform(c);
    point.x = c.v[0];
    point.y = c.v[1];
    proj_context_destroy(ctx);

    return true;
}