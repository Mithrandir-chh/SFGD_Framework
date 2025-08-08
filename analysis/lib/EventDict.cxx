// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME EventDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "classes/Event.hh"
#include "classes/Hit.hh"
#include "classes/ND280SFGDEvent.hh"
#include "classes/ND280SFGDHit.hh"
#include "classes/ND280SFGDTrack.hh"
#include "classes/ND280SFGDVoxel.hh"
#include "classes/ND280SFGDVoxelSet.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_Hit(void *p = nullptr);
   static void *newArray_Hit(Long_t size, void *p);
   static void delete_Hit(void *p);
   static void deleteArray_Hit(void *p);
   static void destruct_Hit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Hit*)
   {
      ::Hit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Hit >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Hit", ::Hit::Class_Version(), "classes/Hit.hh", 8,
                  typeid(::Hit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Hit::Dictionary, isa_proxy, 4,
                  sizeof(::Hit) );
      instance.SetNew(&new_Hit);
      instance.SetNewArray(&newArray_Hit);
      instance.SetDelete(&delete_Hit);
      instance.SetDeleteArray(&deleteArray_Hit);
      instance.SetDestructor(&destruct_Hit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Hit*)
   {
      return GenerateInitInstanceLocal(static_cast<::Hit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Hit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_Event(void *p = nullptr);
   static void *newArray_Event(Long_t size, void *p);
   static void delete_Event(void *p);
   static void deleteArray_Event(void *p);
   static void destruct_Event(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Event*)
   {
      ::Event *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Event >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Event", ::Event::Class_Version(), "classes/Event.hh", 8,
                  typeid(::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Event::Dictionary, isa_proxy, 4,
                  sizeof(::Event) );
      instance.SetNew(&new_Event);
      instance.SetNewArray(&newArray_Event);
      instance.SetDelete(&delete_Event);
      instance.SetDeleteArray(&deleteArray_Event);
      instance.SetDestructor(&destruct_Event);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Event*)
   {
      return GenerateInitInstanceLocal(static_cast<::Event*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Event*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ND280SFGDHit(void *p = nullptr);
   static void *newArray_ND280SFGDHit(Long_t size, void *p);
   static void delete_ND280SFGDHit(void *p);
   static void deleteArray_ND280SFGDHit(void *p);
   static void destruct_ND280SFGDHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ND280SFGDHit*)
   {
      ::ND280SFGDHit *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ND280SFGDHit >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ND280SFGDHit", ::ND280SFGDHit::Class_Version(), "classes/ND280SFGDHit.hh", 12,
                  typeid(::ND280SFGDHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ND280SFGDHit::Dictionary, isa_proxy, 4,
                  sizeof(::ND280SFGDHit) );
      instance.SetNew(&new_ND280SFGDHit);
      instance.SetNewArray(&newArray_ND280SFGDHit);
      instance.SetDelete(&delete_ND280SFGDHit);
      instance.SetDeleteArray(&deleteArray_ND280SFGDHit);
      instance.SetDestructor(&destruct_ND280SFGDHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ND280SFGDHit*)
   {
      return GenerateInitInstanceLocal(static_cast<::ND280SFGDHit*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ND280SFGDHit*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ND280SFGDVoxel(void *p = nullptr);
   static void *newArray_ND280SFGDVoxel(Long_t size, void *p);
   static void delete_ND280SFGDVoxel(void *p);
   static void deleteArray_ND280SFGDVoxel(void *p);
   static void destruct_ND280SFGDVoxel(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ND280SFGDVoxel*)
   {
      ::ND280SFGDVoxel *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ND280SFGDVoxel >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ND280SFGDVoxel", ::ND280SFGDVoxel::Class_Version(), "classes/ND280SFGDVoxel.hh", 11,
                  typeid(::ND280SFGDVoxel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ND280SFGDVoxel::Dictionary, isa_proxy, 4,
                  sizeof(::ND280SFGDVoxel) );
      instance.SetNew(&new_ND280SFGDVoxel);
      instance.SetNewArray(&newArray_ND280SFGDVoxel);
      instance.SetDelete(&delete_ND280SFGDVoxel);
      instance.SetDeleteArray(&deleteArray_ND280SFGDVoxel);
      instance.SetDestructor(&destruct_ND280SFGDVoxel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ND280SFGDVoxel*)
   {
      return GenerateInitInstanceLocal(static_cast<::ND280SFGDVoxel*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ND280SFGDVoxel*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ND280SFGDVoxelSet(void *p = nullptr);
   static void *newArray_ND280SFGDVoxelSet(Long_t size, void *p);
   static void delete_ND280SFGDVoxelSet(void *p);
   static void deleteArray_ND280SFGDVoxelSet(void *p);
   static void destruct_ND280SFGDVoxelSet(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ND280SFGDVoxelSet*)
   {
      ::ND280SFGDVoxelSet *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ND280SFGDVoxelSet >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ND280SFGDVoxelSet", ::ND280SFGDVoxelSet::Class_Version(), "classes/ND280SFGDVoxelSet.hh", 19,
                  typeid(::ND280SFGDVoxelSet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ND280SFGDVoxelSet::Dictionary, isa_proxy, 4,
                  sizeof(::ND280SFGDVoxelSet) );
      instance.SetNew(&new_ND280SFGDVoxelSet);
      instance.SetNewArray(&newArray_ND280SFGDVoxelSet);
      instance.SetDelete(&delete_ND280SFGDVoxelSet);
      instance.SetDeleteArray(&deleteArray_ND280SFGDVoxelSet);
      instance.SetDestructor(&destruct_ND280SFGDVoxelSet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ND280SFGDVoxelSet*)
   {
      return GenerateInitInstanceLocal(static_cast<::ND280SFGDVoxelSet*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ND280SFGDVoxelSet*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ND280SFGDTrack(void *p = nullptr);
   static void *newArray_ND280SFGDTrack(Long_t size, void *p);
   static void delete_ND280SFGDTrack(void *p);
   static void deleteArray_ND280SFGDTrack(void *p);
   static void destruct_ND280SFGDTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ND280SFGDTrack*)
   {
      ::ND280SFGDTrack *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ND280SFGDTrack >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ND280SFGDTrack", ::ND280SFGDTrack::Class_Version(), "classes/ND280SFGDTrack.hh", 11,
                  typeid(::ND280SFGDTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ND280SFGDTrack::Dictionary, isa_proxy, 4,
                  sizeof(::ND280SFGDTrack) );
      instance.SetNew(&new_ND280SFGDTrack);
      instance.SetNewArray(&newArray_ND280SFGDTrack);
      instance.SetDelete(&delete_ND280SFGDTrack);
      instance.SetDeleteArray(&deleteArray_ND280SFGDTrack);
      instance.SetDestructor(&destruct_ND280SFGDTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ND280SFGDTrack*)
   {
      return GenerateInitInstanceLocal(static_cast<::ND280SFGDTrack*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ND280SFGDTrack*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ND280SFGDEvent(void *p = nullptr);
   static void *newArray_ND280SFGDEvent(Long_t size, void *p);
   static void delete_ND280SFGDEvent(void *p);
   static void deleteArray_ND280SFGDEvent(void *p);
   static void destruct_ND280SFGDEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ND280SFGDEvent*)
   {
      ::ND280SFGDEvent *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ND280SFGDEvent >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("ND280SFGDEvent", ::ND280SFGDEvent::Class_Version(), "classes/ND280SFGDEvent.hh", 100,
                  typeid(::ND280SFGDEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ND280SFGDEvent::Dictionary, isa_proxy, 4,
                  sizeof(::ND280SFGDEvent) );
      instance.SetNew(&new_ND280SFGDEvent);
      instance.SetNewArray(&newArray_ND280SFGDEvent);
      instance.SetDelete(&delete_ND280SFGDEvent);
      instance.SetDeleteArray(&deleteArray_ND280SFGDEvent);
      instance.SetDestructor(&destruct_ND280SFGDEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ND280SFGDEvent*)
   {
      return GenerateInitInstanceLocal(static_cast<::ND280SFGDEvent*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::ND280SFGDEvent*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Hit::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Hit::Class_Name()
{
   return "Hit";
}

//______________________________________________________________________________
const char *Hit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hit*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Hit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hit*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Hit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hit*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Hit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hit*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Event*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ND280SFGDHit::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ND280SFGDHit::Class_Name()
{
   return "ND280SFGDHit";
}

//______________________________________________________________________________
const char *ND280SFGDHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDHit*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ND280SFGDHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDHit*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ND280SFGDHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDHit*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ND280SFGDHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDHit*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ND280SFGDVoxel::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ND280SFGDVoxel::Class_Name()
{
   return "ND280SFGDVoxel";
}

//______________________________________________________________________________
const char *ND280SFGDVoxel::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxel*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ND280SFGDVoxel::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxel*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ND280SFGDVoxel::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxel*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ND280SFGDVoxel::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxel*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ND280SFGDVoxelSet::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ND280SFGDVoxelSet::Class_Name()
{
   return "ND280SFGDVoxelSet";
}

//______________________________________________________________________________
const char *ND280SFGDVoxelSet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxelSet*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ND280SFGDVoxelSet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxelSet*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ND280SFGDVoxelSet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxelSet*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ND280SFGDVoxelSet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDVoxelSet*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ND280SFGDTrack::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ND280SFGDTrack::Class_Name()
{
   return "ND280SFGDTrack";
}

//______________________________________________________________________________
const char *ND280SFGDTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDTrack*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ND280SFGDTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDTrack*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ND280SFGDTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDTrack*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ND280SFGDTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDTrack*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ND280SFGDEvent::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *ND280SFGDEvent::Class_Name()
{
   return "ND280SFGDEvent";
}

//______________________________________________________________________________
const char *ND280SFGDEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDEvent*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int ND280SFGDEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDEvent*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ND280SFGDEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDEvent*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ND280SFGDEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ND280SFGDEvent*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Hit::Streamer(TBuffer &R__b)
{
   // Stream an object of class Hit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Hit::Class(),this);
   } else {
      R__b.WriteClassBuffer(Hit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Hit(void *p) {
      return  p ? new(p) ::Hit : new ::Hit;
   }
   static void *newArray_Hit(Long_t nElements, void *p) {
      return p ? new(p) ::Hit[nElements] : new ::Hit[nElements];
   }
   // Wrapper around operator delete
   static void delete_Hit(void *p) {
      delete (static_cast<::Hit*>(p));
   }
   static void deleteArray_Hit(void *p) {
      delete [] (static_cast<::Hit*>(p));
   }
   static void destruct_Hit(void *p) {
      typedef ::Hit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Hit

//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(Event::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Event(void *p) {
      return  p ? new(p) ::Event : new ::Event;
   }
   static void *newArray_Event(Long_t nElements, void *p) {
      return p ? new(p) ::Event[nElements] : new ::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_Event(void *p) {
      delete (static_cast<::Event*>(p));
   }
   static void deleteArray_Event(void *p) {
      delete [] (static_cast<::Event*>(p));
   }
   static void destruct_Event(void *p) {
      typedef ::Event current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Event

//______________________________________________________________________________
void ND280SFGDHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class ND280SFGDHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ND280SFGDHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(ND280SFGDHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ND280SFGDHit(void *p) {
      return  p ? new(p) ::ND280SFGDHit : new ::ND280SFGDHit;
   }
   static void *newArray_ND280SFGDHit(Long_t nElements, void *p) {
      return p ? new(p) ::ND280SFGDHit[nElements] : new ::ND280SFGDHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_ND280SFGDHit(void *p) {
      delete (static_cast<::ND280SFGDHit*>(p));
   }
   static void deleteArray_ND280SFGDHit(void *p) {
      delete [] (static_cast<::ND280SFGDHit*>(p));
   }
   static void destruct_ND280SFGDHit(void *p) {
      typedef ::ND280SFGDHit current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ND280SFGDHit

//______________________________________________________________________________
void ND280SFGDVoxel::Streamer(TBuffer &R__b)
{
   // Stream an object of class ND280SFGDVoxel.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ND280SFGDVoxel::Class(),this);
   } else {
      R__b.WriteClassBuffer(ND280SFGDVoxel::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ND280SFGDVoxel(void *p) {
      return  p ? new(p) ::ND280SFGDVoxel : new ::ND280SFGDVoxel;
   }
   static void *newArray_ND280SFGDVoxel(Long_t nElements, void *p) {
      return p ? new(p) ::ND280SFGDVoxel[nElements] : new ::ND280SFGDVoxel[nElements];
   }
   // Wrapper around operator delete
   static void delete_ND280SFGDVoxel(void *p) {
      delete (static_cast<::ND280SFGDVoxel*>(p));
   }
   static void deleteArray_ND280SFGDVoxel(void *p) {
      delete [] (static_cast<::ND280SFGDVoxel*>(p));
   }
   static void destruct_ND280SFGDVoxel(void *p) {
      typedef ::ND280SFGDVoxel current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ND280SFGDVoxel

//______________________________________________________________________________
void ND280SFGDVoxelSet::Streamer(TBuffer &R__b)
{
   // Stream an object of class ND280SFGDVoxelSet.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ND280SFGDVoxelSet::Class(),this);
   } else {
      R__b.WriteClassBuffer(ND280SFGDVoxelSet::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ND280SFGDVoxelSet(void *p) {
      return  p ? new(p) ::ND280SFGDVoxelSet : new ::ND280SFGDVoxelSet;
   }
   static void *newArray_ND280SFGDVoxelSet(Long_t nElements, void *p) {
      return p ? new(p) ::ND280SFGDVoxelSet[nElements] : new ::ND280SFGDVoxelSet[nElements];
   }
   // Wrapper around operator delete
   static void delete_ND280SFGDVoxelSet(void *p) {
      delete (static_cast<::ND280SFGDVoxelSet*>(p));
   }
   static void deleteArray_ND280SFGDVoxelSet(void *p) {
      delete [] (static_cast<::ND280SFGDVoxelSet*>(p));
   }
   static void destruct_ND280SFGDVoxelSet(void *p) {
      typedef ::ND280SFGDVoxelSet current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ND280SFGDVoxelSet

//______________________________________________________________________________
void ND280SFGDTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class ND280SFGDTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ND280SFGDTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(ND280SFGDTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ND280SFGDTrack(void *p) {
      return  p ? new(p) ::ND280SFGDTrack : new ::ND280SFGDTrack;
   }
   static void *newArray_ND280SFGDTrack(Long_t nElements, void *p) {
      return p ? new(p) ::ND280SFGDTrack[nElements] : new ::ND280SFGDTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_ND280SFGDTrack(void *p) {
      delete (static_cast<::ND280SFGDTrack*>(p));
   }
   static void deleteArray_ND280SFGDTrack(void *p) {
      delete [] (static_cast<::ND280SFGDTrack*>(p));
   }
   static void destruct_ND280SFGDTrack(void *p) {
      typedef ::ND280SFGDTrack current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ND280SFGDTrack

//______________________________________________________________________________
void ND280SFGDEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class ND280SFGDEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(ND280SFGDEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(ND280SFGDEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ND280SFGDEvent(void *p) {
      return  p ? new(p) ::ND280SFGDEvent : new ::ND280SFGDEvent;
   }
   static void *newArray_ND280SFGDEvent(Long_t nElements, void *p) {
      return p ? new(p) ::ND280SFGDEvent[nElements] : new ::ND280SFGDEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_ND280SFGDEvent(void *p) {
      delete (static_cast<::ND280SFGDEvent*>(p));
   }
   static void deleteArray_ND280SFGDEvent(void *p) {
      delete [] (static_cast<::ND280SFGDEvent*>(p));
   }
   static void destruct_ND280SFGDEvent(void *p) {
      typedef ::ND280SFGDEvent current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::ND280SFGDEvent

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 428,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr))->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete (static_cast<vector<int>*>(p));
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] (static_cast<vector<int>*>(p));
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 428,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEND280SFGDVoxelSetgR_Dictionary();
   static void vectorlEND280SFGDVoxelSetgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDVoxelSetgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDVoxelSetgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDVoxelSetgR(void *p);
   static void deleteArray_vectorlEND280SFGDVoxelSetgR(void *p);
   static void destruct_vectorlEND280SFGDVoxelSetgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDVoxelSet>*)
   {
      vector<ND280SFGDVoxelSet> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDVoxelSet>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDVoxelSet>", -2, "vector", 428,
                  typeid(vector<ND280SFGDVoxelSet>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDVoxelSetgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ND280SFGDVoxelSet>) );
      instance.SetNew(&new_vectorlEND280SFGDVoxelSetgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDVoxelSetgR);
      instance.SetDelete(&delete_vectorlEND280SFGDVoxelSetgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDVoxelSetgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDVoxelSetgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDVoxelSet> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDVoxelSet>","std::vector<ND280SFGDVoxelSet, std::allocator<ND280SFGDVoxelSet> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDVoxelSet>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDVoxelSetgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDVoxelSet>*>(nullptr))->GetClass();
      vectorlEND280SFGDVoxelSetgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDVoxelSetgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDVoxelSetgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDVoxelSet> : new vector<ND280SFGDVoxelSet>;
   }
   static void *newArray_vectorlEND280SFGDVoxelSetgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDVoxelSet>[nElements] : new vector<ND280SFGDVoxelSet>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDVoxelSetgR(void *p) {
      delete (static_cast<vector<ND280SFGDVoxelSet>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDVoxelSetgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDVoxelSet>*>(p));
   }
   static void destruct_vectorlEND280SFGDVoxelSetgR(void *p) {
      typedef vector<ND280SFGDVoxelSet> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDVoxelSet>

namespace ROOT {
   static TClass *vectorlEND280SFGDVoxelgR_Dictionary();
   static void vectorlEND280SFGDVoxelgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDVoxelgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDVoxelgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDVoxelgR(void *p);
   static void deleteArray_vectorlEND280SFGDVoxelgR(void *p);
   static void destruct_vectorlEND280SFGDVoxelgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDVoxel>*)
   {
      vector<ND280SFGDVoxel> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDVoxel>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDVoxel>", -2, "vector", 428,
                  typeid(vector<ND280SFGDVoxel>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDVoxelgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ND280SFGDVoxel>) );
      instance.SetNew(&new_vectorlEND280SFGDVoxelgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDVoxelgR);
      instance.SetDelete(&delete_vectorlEND280SFGDVoxelgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDVoxelgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDVoxelgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDVoxel> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDVoxel>","std::vector<ND280SFGDVoxel, std::allocator<ND280SFGDVoxel> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDVoxel>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDVoxelgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDVoxel>*>(nullptr))->GetClass();
      vectorlEND280SFGDVoxelgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDVoxelgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDVoxelgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDVoxel> : new vector<ND280SFGDVoxel>;
   }
   static void *newArray_vectorlEND280SFGDVoxelgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDVoxel>[nElements] : new vector<ND280SFGDVoxel>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDVoxelgR(void *p) {
      delete (static_cast<vector<ND280SFGDVoxel>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDVoxelgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDVoxel>*>(p));
   }
   static void destruct_vectorlEND280SFGDVoxelgR(void *p) {
      typedef vector<ND280SFGDVoxel> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDVoxel>

namespace ROOT {
   static TClass *vectorlEND280SFGDVoxelmUgR_Dictionary();
   static void vectorlEND280SFGDVoxelmUgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDVoxelmUgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDVoxelmUgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDVoxelmUgR(void *p);
   static void deleteArray_vectorlEND280SFGDVoxelmUgR(void *p);
   static void destruct_vectorlEND280SFGDVoxelmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDVoxel*>*)
   {
      vector<ND280SFGDVoxel*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDVoxel*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDVoxel*>", -2, "vector", 428,
                  typeid(vector<ND280SFGDVoxel*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDVoxelmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ND280SFGDVoxel*>) );
      instance.SetNew(&new_vectorlEND280SFGDVoxelmUgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDVoxelmUgR);
      instance.SetDelete(&delete_vectorlEND280SFGDVoxelmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDVoxelmUgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDVoxelmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDVoxel*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDVoxel*>","std::vector<ND280SFGDVoxel*, std::allocator<ND280SFGDVoxel*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDVoxel*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDVoxelmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDVoxel*>*>(nullptr))->GetClass();
      vectorlEND280SFGDVoxelmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDVoxelmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDVoxelmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDVoxel*> : new vector<ND280SFGDVoxel*>;
   }
   static void *newArray_vectorlEND280SFGDVoxelmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDVoxel*>[nElements] : new vector<ND280SFGDVoxel*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDVoxelmUgR(void *p) {
      delete (static_cast<vector<ND280SFGDVoxel*>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDVoxelmUgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDVoxel*>*>(p));
   }
   static void destruct_vectorlEND280SFGDVoxelmUgR(void *p) {
      typedef vector<ND280SFGDVoxel*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDVoxel*>

namespace ROOT {
   static TClass *vectorlEND280SFGDTrackgR_Dictionary();
   static void vectorlEND280SFGDTrackgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDTrackgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDTrackgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDTrackgR(void *p);
   static void deleteArray_vectorlEND280SFGDTrackgR(void *p);
   static void destruct_vectorlEND280SFGDTrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDTrack>*)
   {
      vector<ND280SFGDTrack> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDTrack>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDTrack>", -2, "vector", 428,
                  typeid(vector<ND280SFGDTrack>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDTrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ND280SFGDTrack>) );
      instance.SetNew(&new_vectorlEND280SFGDTrackgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDTrackgR);
      instance.SetDelete(&delete_vectorlEND280SFGDTrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDTrackgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDTrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDTrack> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDTrack>","std::vector<ND280SFGDTrack, std::allocator<ND280SFGDTrack> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDTrack>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDTrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDTrack>*>(nullptr))->GetClass();
      vectorlEND280SFGDTrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDTrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDTrackgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDTrack> : new vector<ND280SFGDTrack>;
   }
   static void *newArray_vectorlEND280SFGDTrackgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDTrack>[nElements] : new vector<ND280SFGDTrack>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDTrackgR(void *p) {
      delete (static_cast<vector<ND280SFGDTrack>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDTrackgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDTrack>*>(p));
   }
   static void destruct_vectorlEND280SFGDTrackgR(void *p) {
      typedef vector<ND280SFGDTrack> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDTrack>

namespace ROOT {
   static TClass *vectorlEND280SFGDTrackmUgR_Dictionary();
   static void vectorlEND280SFGDTrackmUgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDTrackmUgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDTrackmUgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDTrackmUgR(void *p);
   static void deleteArray_vectorlEND280SFGDTrackmUgR(void *p);
   static void destruct_vectorlEND280SFGDTrackmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDTrack*>*)
   {
      vector<ND280SFGDTrack*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDTrack*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDTrack*>", -2, "vector", 428,
                  typeid(vector<ND280SFGDTrack*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDTrackmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ND280SFGDTrack*>) );
      instance.SetNew(&new_vectorlEND280SFGDTrackmUgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDTrackmUgR);
      instance.SetDelete(&delete_vectorlEND280SFGDTrackmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDTrackmUgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDTrackmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDTrack*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDTrack*>","std::vector<ND280SFGDTrack*, std::allocator<ND280SFGDTrack*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDTrack*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDTrackmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDTrack*>*>(nullptr))->GetClass();
      vectorlEND280SFGDTrackmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDTrackmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDTrackmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDTrack*> : new vector<ND280SFGDTrack*>;
   }
   static void *newArray_vectorlEND280SFGDTrackmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDTrack*>[nElements] : new vector<ND280SFGDTrack*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDTrackmUgR(void *p) {
      delete (static_cast<vector<ND280SFGDTrack*>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDTrackmUgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDTrack*>*>(p));
   }
   static void destruct_vectorlEND280SFGDTrackmUgR(void *p) {
      typedef vector<ND280SFGDTrack*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDTrack*>

namespace ROOT {
   static TClass *vectorlEND280SFGDHitgR_Dictionary();
   static void vectorlEND280SFGDHitgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDHitgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDHitgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDHitgR(void *p);
   static void deleteArray_vectorlEND280SFGDHitgR(void *p);
   static void destruct_vectorlEND280SFGDHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDHit>*)
   {
      vector<ND280SFGDHit> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDHit>", -2, "vector", 428,
                  typeid(vector<ND280SFGDHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ND280SFGDHit>) );
      instance.SetNew(&new_vectorlEND280SFGDHitgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDHitgR);
      instance.SetDelete(&delete_vectorlEND280SFGDHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDHitgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDHit> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDHit>","std::vector<ND280SFGDHit, std::allocator<ND280SFGDHit> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDHit>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDHit>*>(nullptr))->GetClass();
      vectorlEND280SFGDHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDHitgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDHit> : new vector<ND280SFGDHit>;
   }
   static void *newArray_vectorlEND280SFGDHitgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDHit>[nElements] : new vector<ND280SFGDHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDHitgR(void *p) {
      delete (static_cast<vector<ND280SFGDHit>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDHitgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDHit>*>(p));
   }
   static void destruct_vectorlEND280SFGDHitgR(void *p) {
      typedef vector<ND280SFGDHit> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDHit>

namespace ROOT {
   static TClass *vectorlEND280SFGDHitmUgR_Dictionary();
   static void vectorlEND280SFGDHitmUgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDHitmUgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDHitmUgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDHitmUgR(void *p);
   static void deleteArray_vectorlEND280SFGDHitmUgR(void *p);
   static void destruct_vectorlEND280SFGDHitmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDHit*>*)
   {
      vector<ND280SFGDHit*> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDHit*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDHit*>", -2, "vector", 428,
                  typeid(vector<ND280SFGDHit*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDHitmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<ND280SFGDHit*>) );
      instance.SetNew(&new_vectorlEND280SFGDHitmUgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDHitmUgR);
      instance.SetDelete(&delete_vectorlEND280SFGDHitmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDHitmUgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDHitmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDHit*> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDHit*>","std::vector<ND280SFGDHit*, std::allocator<ND280SFGDHit*> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDHit*>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDHitmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDHit*>*>(nullptr))->GetClass();
      vectorlEND280SFGDHitmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDHitmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDHitmUgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDHit*> : new vector<ND280SFGDHit*>;
   }
   static void *newArray_vectorlEND280SFGDHitmUgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDHit*>[nElements] : new vector<ND280SFGDHit*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDHitmUgR(void *p) {
      delete (static_cast<vector<ND280SFGDHit*>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDHitmUgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDHit*>*>(p));
   }
   static void destruct_vectorlEND280SFGDHitmUgR(void *p) {
      typedef vector<ND280SFGDHit*> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDHit*>

namespace ROOT {
   static TClass *vectorlEND280SFGDEventgR_Dictionary();
   static void vectorlEND280SFGDEventgR_TClassManip(TClass*);
   static void *new_vectorlEND280SFGDEventgR(void *p = nullptr);
   static void *newArray_vectorlEND280SFGDEventgR(Long_t size, void *p);
   static void delete_vectorlEND280SFGDEventgR(void *p);
   static void deleteArray_vectorlEND280SFGDEventgR(void *p);
   static void destruct_vectorlEND280SFGDEventgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<ND280SFGDEvent>*)
   {
      vector<ND280SFGDEvent> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<ND280SFGDEvent>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<ND280SFGDEvent>", -2, "vector", 428,
                  typeid(vector<ND280SFGDEvent>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEND280SFGDEventgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<ND280SFGDEvent>) );
      instance.SetNew(&new_vectorlEND280SFGDEventgR);
      instance.SetNewArray(&newArray_vectorlEND280SFGDEventgR);
      instance.SetDelete(&delete_vectorlEND280SFGDEventgR);
      instance.SetDeleteArray(&deleteArray_vectorlEND280SFGDEventgR);
      instance.SetDestructor(&destruct_vectorlEND280SFGDEventgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<ND280SFGDEvent> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<ND280SFGDEvent>","std::vector<ND280SFGDEvent, std::allocator<ND280SFGDEvent> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDEvent>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEND280SFGDEventgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<ND280SFGDEvent>*>(nullptr))->GetClass();
      vectorlEND280SFGDEventgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEND280SFGDEventgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEND280SFGDEventgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDEvent> : new vector<ND280SFGDEvent>;
   }
   static void *newArray_vectorlEND280SFGDEventgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<ND280SFGDEvent>[nElements] : new vector<ND280SFGDEvent>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEND280SFGDEventgR(void *p) {
      delete (static_cast<vector<ND280SFGDEvent>*>(p));
   }
   static void deleteArray_vectorlEND280SFGDEventgR(void *p) {
      delete [] (static_cast<vector<ND280SFGDEvent>*>(p));
   }
   static void destruct_vectorlEND280SFGDEventgR(void *p) {
      typedef vector<ND280SFGDEvent> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<ND280SFGDEvent>

namespace {
  void TriggerDictionaryInitialization_libEventDict_Impl() {
    static const char* headers[] = {
"classes/Event.hh",
"classes/Hit.hh",
"classes/ND280SFGDEvent.hh",
"classes/ND280SFGDHit.hh",
"classes/ND280SFGDTrack.hh",
"classes/ND280SFGDVoxel.hh",
"classes/ND280SFGDVoxelSet.hh",
nullptr
    };
    static const char* includePaths[] = {
"/media/disk_b/standard_software/sfgd_framework/analysis/.",
"/media/disk_b/standard_software/sfgd_framework/analysis/src",
"/usr/local/include",
"/usr/local/include/",
"/media/disk_b/standard_software/sfgd_framework/analysis/build/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$classes/ND280SFGDEvent.hh")))  ND280SFGDTrack;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$classes/ND280SFGDEvent.hh")))  ND280SFGDEvent;
class __attribute__((annotate("$clingAutoload$classes/ND280SFGDEvent.hh")))  ND280SFGDHit;
class __attribute__((annotate("$clingAutoload$classes/ND280SFGDEvent.hh")))  ND280SFGDVoxelSet;
class __attribute__((annotate("$clingAutoload$classes/ND280SFGDEvent.hh")))  ND280SFGDVoxel;
class __attribute__((annotate("$clingAutoload$classes/Event.hh")))  Hit;
class __attribute__((annotate("$clingAutoload$classes/Event.hh")))  Event;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "classes/Event.hh"
#include "classes/Hit.hh"
#include "classes/ND280SFGDEvent.hh"
#include "classes/ND280SFGDHit.hh"
#include "classes/ND280SFGDTrack.hh"
#include "classes/ND280SFGDVoxel.hh"
#include "classes/ND280SFGDVoxelSet.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"Event", payloadCode, "@",
"Hit", payloadCode, "@",
"ND280SFGDEvent", payloadCode, "@",
"ND280SFGDHit", payloadCode, "@",
"ND280SFGDTrack", payloadCode, "@",
"ND280SFGDVoxel", payloadCode, "@",
"ND280SFGDVoxelSet", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libEventDict() {
  TriggerDictionaryInitialization_libEventDict_Impl();
}
