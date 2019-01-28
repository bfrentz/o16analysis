#ifndef CONFIG_H
#define CONFIG_H

///Indicates the evt files come from the NSCL "Ring" Buffer.
#define RING_BUFFER false 
///Indicates the evt files come from a VM-USB crate.
#define VM_USB true

//Uncomment the following line to use a manually specified buffer size.
// Most setups use the default buffer size and this line can be ignored.
//#define BUFFER_SIZE 4096 ///< Manually defined buffer size.

///The list of modules to be unpacked.
/**Modules must be listed in the order that they are to be unpacked. 
 * (Be sure to continue the define line with the '\' character.)
 * If multiple modules of the same type are used a line for each
 * module is required. See modules/include for available modules.
 * The predefine is read by cmake to build a library specific to
 * each implementation and created a vector of pointers to the
 * ReadEvent() methods of each module class.
 * 
 * New modules can be added with the following required prototype:
 * \code 
 * ReadEvent(nsclBuffer*,eventData*,bool)
 * \endcode
 * See modules/include for examples.
 */
#define MODULE_LIST(MODULE) \
	MODULE(Mesytec_MQDC32) \
	MODULE(Caen_General) \
	MODULE(ASIC_XLM)

//MESYTEC MQDC-32 MODULE
#define MESYTEC_MQDC32_ID 5

//MESYTEC MADC MODULE
#define MESYTEC_MADC_ID 5
/* Mesytec resolution vlaue
 * 0 = 2k, 800ns conv.
 * 1 = 4k, 1.6us conv.
 * 2 = 4k, hires, 3.2us conv.
 * 3 = 8k, 1.6us conv.
 * 4 = 8k, hires, 3.2us conv.
 */
#define MESYTEC_RESOLUTION 4

//CAEN TDC MODULE
#define CAEN_GEO 17
#define CAEN_CRATE 0

//ASIC number of chip boards (32ch each)
#define ASIC_NUM_CB 9
#define ASIC_E_SLOT 0
#define ASIC_T_SLOT 1
#define ASIC_MAX_HIT 32

#endif

