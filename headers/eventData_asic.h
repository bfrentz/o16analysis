#ifndef EVENTDATA_H
#define EVENTDATA_H

#include "evt_config.h"
#include "TObject.h"
#include "TRandom3.h"

class eventData : public TObject
{
	private:
		
		// CAEN TDC 32ch
		UShort_t tdc[32]; //!temporary
//		UShort_t tdc[32];

		// Mesytec MQDC-32
		UShort_t qdc[32]; //!temporary
//		UShort_t qdc[32];
		//Float_t sili_en[6];

		// HINP asic chips
//		UShort_t asic_nhit; //!temporary
//		UShort_t asic_ncb[ASIC_MAX_HIT]; //!temporary
//		UShort_t asic_nch[ASIC_MAX_HIT]; //!temporary
//		UShort_t asic_hite[ASIC_MAX_HIT]; //!temporary
//		UShort_t asic_hitt[ASIC_MAX_HIT]; //!temporary
//		UShort_t asic_e[ASIC_NUM_CB][32]; //!temporary
//		UShort_t asic_t[ASIC_NUM_CB][32]; //!temporary
		UShort_t asic_nhit;
		UShort_t asic_ncb[ASIC_MAX_HIT];
		UShort_t asic_nch[ASIC_MAX_HIT];
		UShort_t asic_hite[ASIC_MAX_HIT];
		UShort_t asic_hitt[ASIC_MAX_HIT];
		UShort_t asic_e[ASIC_NUM_CB][32];
		UShort_t asic_t[ASIC_NUM_CB][32];

	public:
		///Default Constructor.
		eventData();
		///Reset the current values.
		void Reset();
		///Record the value from the specified location.
		void SetValue(int crate, int slot, int ch, int value);
		///Record the data array for ASIC
		void SetValues(int crate, int slot, int nhit, int ch[], int evalue[], int tvalue[]);
		///Calibrates the data from one event.
		void Calibrate();



	ClassDef(eventData,1);
};

#endif


